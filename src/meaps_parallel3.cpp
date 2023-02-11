// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "deborder.h"
#include "distribuer.h"
#include "utils_newton_methods.h"

using namespace Rcpp;
using namespace RcppParallel;

inline std::vector<double> repartir_actifs(std::vector<double> dispo, 
                                           std::vector<double> od,
                                           double fuite, 
                                           double actifs) {
  
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
      eps = abs(new_cref - c_ref);
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

struct Meaps_parallel2 : public Worker
{   
  // Les inputs.
  const RMatrix<int> rkdist;
  const RVector<double> emplois;
  const RVector<double> actifs;
  const RMatrix<double> modds;
  const RVector<double> f;
  const RMatrix<int> shuf;
  const std::size_t N;
  const std::size_t K;
  const std::size_t Ns;
  
  // La valeur accumulée.
  std::vector<double> liaisons;
  
  // Constructeur 1.
  Meaps_parallel2(const IntegerMatrix rkdist,
                  const NumericVector emplois,
                  const NumericVector actifs,
                  const NumericMatrix modds,
                  const NumericVector f,
                  const IntegerMatrix shuf,
                  const std::size_t N,
                    const std::size_t K,
                    const std::size_t Ns) : rkdist(rkdist), emplois(emplois), actifs(actifs), 
                  modds(modds), f(f), shuf(shuf), N(N), K(K), Ns(Ns),
                  liaisons() {
    liaisons.resize(N*K, 0.0);
  }
  
  // Constructeur 2 pour le split.
  Meaps_parallel2(const Meaps_parallel2 &meaps, Split) : rkdist(meaps.rkdist), emplois(meaps.emplois), 
      actifs(meaps.actifs), modds(meaps.modds), f(meaps.f), 
      shuf(meaps.shuf), N(meaps.N), K(meaps.K), Ns(meaps.Ns), liaisons() {
    liaisons.resize(N*K, 0.0);
  }
  
  // Accumule les résultats sur la plage souhaitée.
  void operator()(std::size_t begin, std::size_t end) {
    
    // std::size_t N = rkdist.nrow(),
    //   K = rkdist.ncol(),
    //   Ns = shuf.ncol();
    
    Rcout << "begin = " << begin << "\n";
    Rcout << "end = " << end << "\n";
    
    //Initialisation de la matrice de résultat.
    // for (std::size_t i = 0; i < N; ++i){
    //   for (std::size_t j = 0; j < K; ++j) {
    //     liaisons(i,j) = 0.0;
    //   }
    // }
    return;
    for (auto s = begin; s < end; ++s) {
      
      
      
      //RMatrix<int>::Row ishuf = shuf.row(s);
      // for (std::size_t k = 0; k < Ns; ++k) {
      //   ishuf[k] -= 1L;
      // }
      
      std::vector<std::size_t> ishuf(Ns);
      for (std::size_t k = 0; k < Ns; ++k) {
        ishuf[k] = shuf(s, k) - 1L;
      }
      
      
      
      // Le vecteur shuf peut être plus long que le nombre de lignes de rkdist s'il fait repasser plusieurs fois
      // la même ligne d'actifs. Dans ce cas, on compte la fréquence de passage de chaque ligne et l'on divise le
      // poids de la ligne par cette fréquence.
      std::vector<int> freq_actifs(N, 0L);
      for (auto i: ishuf) {
        freq_actifs[i]++;
      }
      
      
      
      std::vector<double> emp(K);
      for (std::size_t k = 0; k < K; ++k) {
        emp[k] = emplois[k];
      }
      
      
      
      for (auto i: ishuf) {
        
        RMatrix<int>::Row rki = rkdist.row(i);
        RMatrix<double>::Row odds = modds.row(i);
        
        int temp;
        std::size_t k_valid = 0;
        std::vector<int> arrangement(K);
        for (std::size_t k = 0; k < K; ++k) {
          temp = rki[k];
          if (R_IsNA(temp) == false) { 
            arrangement[temp - 1L] = k;
            k_valid++;
          }
        }
        arrangement.resize(k_valid);
        
        std::vector<double> dispo (k_valid), 
        od(k_valid), 
        repartition(k_valid);
        for (std::size_t j = 0; j < k_valid; ++j) {
          dispo[j] = emp[arrangement[j]]; 
          od[j] = odds[arrangement[j]]; 
        }
        
        // Choix d'une limite basse pour la fuite.
        double fuite = std::max(1e-3, f[i]);
        double actifs_inzone = (1 - fuite) * actifs[i] / freq_actifs[i];
        
        repartition = repartir_actifs(dispo, od, fuite, actifs_inzone);
        
        // Inscription des résultats locaux dans la perspective globale.
        
        
        for(std::size_t j = 0; j < k_valid ; ++j) {
          emp[arrangement[j]] -= repartition[j];
          liaisons[i * K + arrangement[j]] += repartition[j];
          Rcout << liaisons[i * K + arrangement[j]] << "   ";
        }
      }
    }
    
    // for (std::size_t i = 0; i < N; ++i) {
    //   for (std::size_t j = 0; j < K; ++j) {
    //     liaisons(i,j) = resultat[i][j];
    //   }
    // }
  }
  
  void join(const Meaps_parallel2 &rhs) {
    std::size_t N = rkdist.nrow(),
      K = rkdist.ncol();
    Rcout << "join liaisons(1,1) = " << liaisons[0] << "\n";
    Rcout << "join rhs(1,1) = " << rhs.liaisons[0]<< "\n";
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < K; ++j) {
        liaisons[i * K + j] += rhs.liaisons[i * K + j];
      }
    }
  }
  
};

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
  
  int N = rkdist.nrow(),
    K = rkdist.ncol(),
    Nboot = shuf.nrow();
  
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
  
  NumericMatrix liaisons2(N,K);
  // Attention : calage des emplois sur le nombre d'actifs.
  NumericVector emp = emplois / sum(emplois) * sum(actifs * (1 - f)); 
  // Déclaration de l'objet Meaps_parallel.
  
  Meaps_parallel2 meaps_res(rkdist, emp, actifs, modds, f, shuf, N, K, shuf.ncol());
  
  // Call.
  parallelReduce(0, Nboot, meaps_res, 1L);
  
  //Retour de la matrice liaisons cumulées divisées par Nboot.
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < K; ++j) {
      liaisons2(i,j) = meaps_res.liaisons[i * K + j] / Nboot;
    }
  }
  
  return liaisons2;
}


