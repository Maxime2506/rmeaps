// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <valarray>
#ifdef _OPENMP
# include <omp.h>
#endif
#include "repartir_actifs.h"

using namespace Rcpp;

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
            Nboot = shuf.nrow(),
            NK = N * K;
  
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
  omp_set_num_threads(1L);
  #endif
  
  // Conversion en objet C++.
  std::vector<int> rkdistcpp = as<std::vector<int>>(rkdist);
  std::vector<double> moddscpp = as<std::vector<double>>(modds);
  std::vector<int> ishufcpp = as<std::vector<int>>(shuf);
  std::vector<double> fcpp = as<std::vector<double>>(f);
  std::vector<double> actifscpp = as<std::vector<double>>(actifs);
  
  // Passage aux valarrays.
  std::valarray<int> ranking(rkdistcpp.data(), NK);
  const std::valarray<double> odds(moddscpp.data(), NK);
  std::valarray<int> ishuf(ishufcpp.data(), Ns*Nboot);
  // Décalage des rangs d'une unité.
  ranking -= 1L;
  ishuf -= 1L;
  
  // Initialisation du résultat.
  // std::vector<double> liaisons(N*K);
  double liaisons[NK] = {}; // essai avec un array.
  
  // Conversion de l'emploi en c++ et calage.
  std::vector<double> emploisinitial = as<std::vector<double>>(emplois);
  // Attention : calage des emplois sur le nombre d'actifs.
  double cale = sum(actifs * (1 - f)) / sum(emplois); 
  for (auto& k: emploisinitial) { k *= cale; }
  
  // Lancement du bootstrap.
 
  #ifdef _OPENMP
  #pragma omp parallel for shared(Nboot, N, K, ishuf, emploisinitial, ranking, odds, fcpp, actifscpp) \
   // reduction (+ : liaisons[:NK])
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
    Rcout << "freqactif \n";
    for (int i =0; i<Ns;++i) Rcout << freq_actifs[i] << " / ";
    Rcout << "\n";
      
      
    for (auto i: theshuf) {
      
      std::valarray<int> rki = ranking[std::slice(i, K, N)];
      std::valarray<double> lignodds = odds[std::slice(i, K, N)];
      
      int temp;
      std::size_t k_valid = 0;
      std::vector<int> arrangement(K);
      for (int k = 0; k < K; ++k) {
        temp = rki[k];
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
        Rcout << "k = " << k << " rep = " << repartition[k] << "\n";
        emp[arrangement[k]] -= repartition[k];
        liaisons[i + N * arrangement[k]] += repartition[k]; // Attention : ordre de remplissage en "ligne".
      }
    }
  }
  
  // Passage d'un array à NumericMatrix, divisées par le nombre de tirage Nboot.
  NumericMatrix resultat(N, K);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < K; ++j) {
      resultat(i,j) = liaisons[i + N * j] / Nboot ;
    } 
  }
  return resultat;
}



