#ifdef _OPENMP
# include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <valarray>
#include "repartir_actifs.h"

using namespace Rcpp;

//' MEAPS en calculant la tension sur les opportunités, 
//' c'est-à-dire le rang moyen sur les shufs juste avant saturation.
//' @param rkdist La matrice des rangs dans lequel les colonnes j sont passées en revue pour chacune des lignes i.
//' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes).
//' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
//' @param modds La matrice des odds modifiant la chance d'absorption de chacun des sites j pour des résidents en i.
//' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude.
//' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs.
//' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto.
//' @param progress Ajoute une barre de progression. Default : true.
//' @param normalisation Calage des emplois disponibles sur le nombre d'actifs travaillant sur la zone. Default : false.
//' @param fuite_min Seuil minimal pour la fuite d'un actif. Doit être supérieur à 0. Défault = 1e-3.
//' @param seuil_newton Seuil relatif pour la convergence par newton du calage de la probabilité d'absorption. Défault = 1e-6.
//' @param seuil_dispo seuil absolu à partir duquel une opportunité est réputée saturée. Défaut 0.1
//' 
//' @return renvoie une matrice avec les estimations du nombre de trajets de i vers j.
// [[Rcpp::export]]
List meaps_tension(IntegerMatrix rkdist,
                            NumericVector emplois,
                            NumericVector actifs,
                            NumericMatrix modds,
                            NumericVector f,
                            IntegerMatrix shuf,
                            int nthreads = 0,
                            bool progress = true,
                            bool normalisation = false,
                            double fuite_min = 1e-3,
                            double seuil_newton = 1e-6,
                            double seuil_dispo = 0.1) {
  
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
  int ntr = nthreads;
  if (ntr == 0) {ntr = omp_get_max_threads();}
  if(ntr>omp_get_max_threads()) {ntr = omp_get_max_threads();}
  if(progress==TRUE) REprintf("Nombre de threads = %i\n", ntr);
  //omp_set_num_threads(ntr);
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
  std::vector<double> liaisons((N + 2L)*K, 0.0); // N*K valeurs de flux, puis K nombre de fois un site est complet, puis K rang complétude.
  
  
  // Conversion de l'emploi en c++ et calage.
  std::vector<double> emploisinitial = as<std::vector<double>>(emplois);
  // Attention : calage des emplois sur le nombre d'actifs.
  if (normalisation) {
    double cale = sum(actifs * (1 - f)) / sum(emplois); 
    for (auto& k: emploisinitial) { k *= cale; }
  }
  // Lancement du bootstrap.
  Progress p(Nboot * Ns, progress);
#ifdef _OPENMP
#pragma omp declare reduction(vsum : std::vector<double> :              \
  std::transform(omp_out.begin(), omp_out.end(),                        \
                 omp_in.begin(), omp_out.begin(), std::plus<double>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp parallel for num_threads(ntr)                               \
    shared(Nboot, Ns, N, K, ishuf, emploisinitial, ranking, odds, fcpp, actifscpp) \
    reduction (vsum : liaisons)
#endif
    for (int iboot = 0; iboot < Nboot; ++iboot) {
      
      std::valarray<int> theshuf = ishuf[std::slice(iboot, Ns, Nboot)];
      
      // Initialisation de l'emploi disponible.
      std::vector<double> emp(K);
      std::copy(emploisinitial.begin(), emploisinitial.end(), emp.begin());
      
      // Le vecteur shuf peut être plus long que le nombre de lignes de rkdist s'il fait repasser plusieurs fois
      // la même ligne d'actifs. Dans ce cas, on compte la fréquence de passage de chaque ligne et l'on divise le
      // poids de la ligne par cette fréquence.
      std::vector<int> freq_actifs(N, 0L);
      for (auto i: theshuf) {
        freq_actifs[i]++;
      }
      
      double rang_actifs = 0.0;
      
      for (auto i: theshuf) {
        
        // Increment progress_bar.
        p.increment();
        
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
        double fuite = std::max(fuite_min, fcpp[i]);
        // Nombre d'actifs en emplois dans la zone repartis en freq_actif paquets.
        double actifs_inzone = (1 - fuite) * actifscpp[i] / freq_actifs[i];
        rang_actifs += actifs_inzone;
        
        repartition = repartir_actifs(dispo, od, fuite, actifs_inzone, seuil_newton);
        
        // Inscription des résultats locaux dans la perspective globale.
        for(std::size_t k = 0; k < k_valid ; ++k) {
          if ((dispo[k]-repartition[k]) < seuil_dispo && dispo[k] > seuil_dispo) { 
            liaisons[ N*K + arrangement[k] ] += 1.0; 
            liaisons[ (N + 1L) * K + arrangement[k] ] += rang_actifs;
            }
          emp[ arrangement[k] ] -= repartition[k];
          liaisons[ i + N * arrangement[k] ] += repartition[k]; // Attention : ordre de remplissage en "ligne".
        }
      }
    }
    
  // Passage d'un vecteur à une NumericMatrix, divisées par le nombre de tirage Nboot.
  NumericMatrix resultat(N, K);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < K; ++j) {
      resultat(i,j) = liaisons[i + N * j] / Nboot ;
    } 
  }
  
  double tot_actifs = sum(actifs);
  NumericVector tension(K, 0.0);
  for (int j = 0; j < K; ++j) {
    tension[j] = ( (Nboot - liaisons[ N * K + j ]) * tot_actifs + liaisons[ (N + 1L) * K + j] ) / Nboot ;
  }
  
  return List::create(_["flux"]= resultat, 
                      _["tension"] = tension);
}


