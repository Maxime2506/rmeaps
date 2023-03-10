// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <valarray>
#ifdef _OPENMP
# include <omp.h>
#endif
#include "deborder.h"
#include "distribuer.h"
#include "utils_newton_methods.h"

using namespace Rcpp;

inline std::vector<double> repartir_actifs(std::vector<double>& dispo, 
                                           std::vector<double>& od,
                                           double& fuite, 
                                           double& actifs) {
  
  int k_valid = dispo.size();
  std::vector<double> repartition(k_valid);
  // Calcul de la place disponible total sur la ligne.
  double tot = 0;
  for (auto& d: dispo) {
    tot += d;
  }
  
  if (tot < 1e-2 || k_valid == 1) {
    // Seuil où l'emploi disponible est suffisamment petit pour s'épargner des calculs inutiles et fragiles.
    repartition = deborder(dispo, actifs);
    
  } else {
    // Calcul par la méthode de Newton de la chance d'absorption de référence compatible avec la fuite.
    double p_ref, c_ref, new_cref, eps;
    
    p_ref = 1  - pow(fuite, 1 / tot);
    c_ref = p_ref / (1 - p_ref);// Chance d'absorption de référence. Calcul initial non calé sur la fuite.
    std::vector<double> proportions(k_valid);
    
    int compteur = 0;
    do {
      new_cref = c_ref - (log_fuite(c_ref, dispo, od) + log(fuite))/ d_logfuite(c_ref, dispo, od); 
      eps = std::abs(new_cref - c_ref);
      c_ref = new_cref;
      compteur++;
      if (compteur > 1000) {
        Rcout << "Le calcul par la méthode de Newton de la chance d\'absorption n\'a pas convergé.\n";
        eps = 0; // L'estimation est poursuivie malgré tout...
      }
    } while (eps > 1e-6);
    
    std::vector<double> c_abs(k_valid);
    for (int j = 0; j < k_valid; ++j) {
      c_abs[j] = od[j] * c_ref;
    }
    // Calcul des proba d'arrivées sur chacun des sites (dépend du chemin, != p_abs).
    double logpass = 0.0, 
      logfuit;
    for (int j = 0; j < k_valid; ++j) {
      logfuit = - dispo[j] * log(1 + c_abs[j]); // proba conditionnelle en log de fuir j une fois arrivée jusqu'à j.
      proportions[j] = exp(logpass) * (1 - exp(logfuit));
      logpass += logfuit; // proba en log d'arriver jusqu'à j+1, calculé à partir de la proba d'arriver jusqu'en j.
    }
    
    repartition = distribuer(dispo, proportions, actifs);
    
  }
  return repartition;
} 

//' @param rkdist La matrice des rangs dans lequel les colonnes j sont passées en revue pour chacune des lignes i.
//' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes).
//' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
//' @param modds La matrice des odds modifiant la chance d'absorption de chacun des sites j pour des résidents en i.
//' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude.
//' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs.
//' 
//' @return renvoie une matrice avec les estimations du nombre de trajets de i vers j.
// [[Rcpp::export]]
NumericMatrix meaps_bootstrap2(IntegerMatrix rkdist,
                               NumericVector emplois,
                               NumericVector actifs,
                               NumericMatrix modds,
                               NumericVector f,
                               IntegerMatrix shuf) {
  
  const int N = rkdist.nrow(),
            K = rkdist.ncol(),
            Ns = shuf.ncol(),
            Nboot = shuf.nrow();
  
  // Quelques vérifs préalables.
  if (emplois.size() != K) {
    stop("Le vecteur emplois et le nb de colonnes de la matrice rkdist ne correspondent pas.");
  }
  if (f.size() != N) {
    stop("Le vecteur fuite et le nb de lignes de la matrice rkdist ne correspondent pas.");
  }
  if (actifs.size() != N) {
    stop("Le vecteur actifs et le nb de lignes de la matrice rkdist ne correspondent pas.");
  }
  if (modds.ncol() != K) {
    stop("La matrice modds et la matrice rkdist n\'ont pas le même nombre de colonnes.");
  }
  if (modds.nrow() != N) {
    stop("La matrice modds et la matrice rkdist n\'ont pas le même nombre de lignes.");
  }
  
  #ifdef _OPENMP
  Rcout << "Test : Openmp est bien défini !";
  omp_set_num_threads(1L);
  #endif
  
  // Conversion en objet C++.
  std::vector<int> rkdistcpp = as<std::vector<int>>(rkdist);
  std::vector<double> moddscpp = as<std::vector<double>>(modds);
  std::vector<int> ishufcpp = as<std::vector<int>>(shuf);
  std::vector<double> fcpp = as<std::vector<double>>(f);
  std::vector<double> actifscpp = as<std::vector<double>>(actifs);
  
  // Passage aux valarrays.
  std::valarray<int> ranking(rkdistcpp.data(), N*K);
  const std::valarray<double> odds(moddscpp.data(), N*K);
  std::valarray<int> ishuf(ishufcpp.data(), Ns*Nboot);
  // Décalage des rangs d'une unité.
  ranking -= 1L;
  ishuf -= 1L;
  
  // Initialisation du résultat.
  double liaisons[N*K]; // array semble être requis pour reduction dans openmp.
  
  // Conversion de l'emploi en c++ et calage.
  std::vector<double> emploisinitial = as<std::vector<double>>(emplois);
  // Attention : calage des emplois sur le nombre d'actifs.
  double cale = sum(actifs * (1 - f)) / sum(emplois); 
  for (auto& k: emploisinitial) { k *= cale;}
  
  // Lancement du bootstrap.
  #ifdef _OPENMP
  #pragma omp parallel for shared(Nboot, N, K, ishuf, emploisinitial, ranking, odds, fcpp, actifscpp) \
    reduction (+:liaisons[:N*K])
  #endif
  for (int iboot = 0; iboot < Nboot; ++iboot) {
    
    std::valarray<int> theshuf = ishuf[std::slice(iboot, Ns, Nboot)];
    
    // Initialisation de l'emploi disponible.
    std::vector<double> emp(emploisinitial);
    
    // Le vecteur shuf peut être plus long que le nombre de lignes de rkdist s'il fait repasser plusieurs fois
    // la même ligne d'actifs. Dans ce cas, on compte la fréquence de passage de chaque ligne et l'on divise le
    // poids de la ligne par cette fréquence.
    std::vector<int> freq_actifs(N, 0L);
    for (auto i: theshuf) {
      freq_actifs[i]++;
    }
    
    for (auto i: theshuf) {
      
      std::valarray<int> rki = ranking[std::slice(i, K, N)];
      std::valarray<double> lignodds = odds[std::slice(i, K, N)];
      
      int temp;
      std::size_t k_valid = 0;
      std::vector<int> arrangement(K);
      for (int k = 0; k < K; ++k) {
        temp = ranking[k];
        if (R_IsNA(temp) == false) { 
          arrangement[temp] = k; //Attention : rkdist rank à partir de 1. Décalage déjà pris en compte.
          k_valid++;
        }
      }
      arrangement.resize(k_valid);
      
      std::vector<double> dispo (k_valid), od(k_valid), repartition(k_valid);
      for (std::size_t k = 0; k < k_valid; ++k) {
        dispo[k] = emp[arrangement[k]]; 
        od[k] = lignodds[arrangement[k]]; 
      }
      
      // Choix d'une limite basse pour la fuite.
      double fuite = std::max(1e-3, fcpp[i]);
      // Nombre d'actifs en emplois dans la zone repartis en freq_actif paquets.
      double actifs_inzone = (1 - fuite) * actifscpp[i] / freq_actifs[i];
      
      repartition = repartir_actifs(dispo, od, fuite, actifs_inzone);
      
      // Inscription des résultats locaux dans la perspective globale.
      for(std::size_t k = 0; k < k_valid ; ++k) {
        emp[arrangement[k]] -= repartition[k];
        liaisons[i * N + arrangement[k]] += repartition[k]; // Attention : ordre de remplissage en "ligne".
      }
    }
  }
  
  // Passage d'un array à NumericMatrix, divisées par le nombre de tirage Nboot.
  NumericMatrix resultat(N, K);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < K; ++j) {
      resultat(i,j) = liaisons[i*N + j]/ Nboot ;
    } 
  }
  return resultat;
}



