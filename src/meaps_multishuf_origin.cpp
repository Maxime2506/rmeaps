#ifdef _OPENMP
# include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "repartir_actifs.h"
#include "fcts_penal.h"

using namespace Rcpp;
//' La fonction MEAPS sur plusieurs shufs
 //' version optimisée
 //' @param jr_dist La matrice des rangs dans lequel les colonnes j sont passées en revue pour chacune des lignes i.
 //' @param p_dist La matrice des rangs dans lequel les colonnes j sont passées en revue pour chacune des lignes i.
 //' @param xr_dist La matrice des rangs dans lequel les colonnes j sont passées en revue pour chacune des lignes i.
 //' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes).
 //' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
 //' @param xr_odds La matrice des odds modifiant la chance d'absorption de chacun des sites j pour des résidents en i.
 //' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude.
 //' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs.
 //' @param attraction l'attraction
 //' @param parametres les paramètres
 //' @param mode Choix du rôle de l'emploi disponible au cours du processus. Default : continu. Autre choix : discret, subjectif_c ou _d...
 //' @param odds_subjectifs Attractivité des sites proches des actifs, pour le mode defini. default : null.
 //' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto.
 //' @param progress Ajoute une barre de progression. Default : true.
 //' 
 //' @return renvoie une matrice avec les estimations du nombre de trajets de i vers j.
 // [[Rcpp::export]]
 List multishuf_origin_cpp(const IntegerVector jr_dist, 
                          const IntegerVector p_dist, 
                          const NumericVector xr_dist, 
                          const NumericVector emplois,
                          const NumericVector actifs, 
                          const NumericVector fuites, 
                          const IntegerMatrix shuf,
                          const std::string attraction,
                          const NumericVector parametres,
                          const NumericVector xr_odds,
                          const std::string mode = "continu",
                          const Nullable<NumericVector> oddssubjectifs = R_NilValue,
                          const int nthreads = 0, 
                          const bool verbose = true) {
   
   const double SEUIL_NEWTON = 1e-6;
   
   const int N = actifs.size(), K = emplois.size(), L = xr_dist.size(), Nx  = xr_dist.size(),
     Ns = shuf.ncol(), Nboot = shuf.nrow();
   
   // Passage explicite en std::vector pour rendre les vecteurs thread safe (ts_)(nécessaire pour openmp dans la macro).
   std::vector<double> ts_emplois = as<std::vector<double>>(emplois);
   std::vector<double> ts_fuite = as<std::vector<double>>(fuites);
   const std::vector<double> ts_actifs = as<std::vector<double>>(actifs);
   
   const std::vector<double> ts_parametres = as< std::vector<double> >(parametres);
   
   const std::vector<int> ts_jr_dist = as< std::vector<int> >(jr_dist);
   const std::vector<int> ts_p_dist = as< std::vector<int> >(p_dist);
   const std::vector<double> ts_xr_dist = as< std::vector<double> >(xr_dist);
   
   std::vector<double> ts_xr_odds(L, 1);
   
   if (attraction == "marche") {
     for (std::size_t k = 0; k < L; ++k) {
       ts_xr_odds[k]  = marche(ts_xr_dist[k], parametres[0], parametres[1]);
     }}
   
   if (attraction == "marche_liss") {
     for (std::size_t k = 0; k < L; ++k) {
       ts_xr_odds[k]  = marche_liss(ts_xr_dist[k], parametres[0], parametres[1]);
     }}
   
   if (attraction == "double_marche_liss") {
     for (std::size_t k = 0; k < L; ++k) {
       ts_xr_odds[k]  = double_marche_liss(ts_xr_dist[k], parametres[0], parametres[1], parametres[2], parametres[3]);
     }}
   
   if (attraction == "decay") {
     for (std::size_t k = 0; k < L; ++k) {
       ts_xr_odds[k]  = decay(ts_xr_dist[k], parametres[0], parametres[1]);
     }}
   
   if (attraction == "logistique") {
     for (std::size_t k = 0; k < L; ++k) {
       ts_xr_odds[k] = logistique(ts_xr_dist[k], parametres[0], parametres[1], parametres[2]);
     }}
   
   if (attraction == "odds") {
     for (std::size_t k = 0; k < L ; ++k) {
       ts_xr_odds[k] =  xr_odds[k]; // dans ce cas, on suppose que les odds sont passés comme avant (pas en log)
     } 
   }
   
   std::vector<double> odsub;
   if ((mode == "subjectif_c"|| mode == "subjectif_d") && oddssubjectifs.isNotNull()) {
     odsub = as< std::vector<double> >(oddssubjectifs);
   }
   
   // Les rangs dans shuf commencent à 1.
   std::vector< std::vector<int> > ishuf(Nboot, std::vector<int>(Ns));
   for (std::size_t i = 0; i < Nboot; ++i) {
     for (std::size_t j = 0; j < Ns; ++j) {
       ishuf[i][j] = shuf(i, j);
       if (ishuf[i][j] >= N) {
         stop("Les rangs de shufs vont au-delà de la taille du vecteur actifs.");
       }
     }
   }
   
   // Le vecteur shuf peut être plus long que le nombre de lignes de rkdist s'il fait repasser plusieurs fois
   // la même ligne d'actifs. Dans ce cas, on compte la fréquence de passage de chaque ligne et l'on divise le
   // poids de la ligne par cette fréquence.
   
   std::vector<int> freq_actifs(N, 0L);
   for (auto i : ishuf[0]) {
     freq_actifs[i]++;
   }
   
#ifdef _OPENMP
   int ntr = nthreads;
   if (ntr == 0) {
     ntr = omp_get_max_threads();
   }
   if (ntr > omp_get_max_threads()) {
     ntr = omp_get_max_threads();
   }
   if (verbose == TRUE) REprintf("Nombre de threads = %i\n", ntr);
#endif
   
   // Initialisation du résultat.
   // Lancement du bootstrap.
   Progress p(Nboot*Ns, verbose);
   
   std::size_t NNboot = floor( Nboot / ntr );
   if(NNboot == 0) NNboot = 1;
   // Un vecteur représentant la matrice des flux groupés.
   std::vector< std::vector<float> > liaisons(ntr, std::vector<float> (Nx));
   for(size_t Iboot = 0; Iboot < ntr; ++Iboot) {
     std::fill( liaisons[Iboot].begin(), liaisons[Iboot].end(), 0 );
   }
   
#pragma omp parallel num_threads(ntr) 
{
#pragma omp for schedule(static)
  for(size_t Iboot = 0; Iboot < ntr; ++Iboot) {
    size_t max_i = (Iboot+1)*NNboot;
    if(Iboot==ntr-1) max_i = Nboot;
    size_t min_i = Iboot*NNboot;
    if(min_i<max_i) {
      for (size_t iboot = min_i; iboot < max_i; ++iboot) {
        // Initialisation de l'ordre de départ des actifs et de l'emploi
        // disponible au début.
        
        // Initialisation de l'emploi disponible.
        std::vector<double> emploisdisp(ts_emplois); // deep copy.
        
        for (auto i: ishuf[iboot]) {
          
          // Increment progress_bar.
          p.increment();
          
          auto debut = ts_p_dist[i], fin = ts_p_dist[i +1L];
          auto k_valid = fin - debut;
          
          // Cas particulier où il n'y a aucune destination valide.
          // ATTENTION je crois que les sites sans places libres peuvent ne pas être éliminées pour la proc répartir.
          //if (k_valid == 0) continue;
          
          std::vector<double> placeslibres (k_valid), attractivites(k_valid), repartition(k_valid);
          for (auto k = 0; k < k_valid; ++k) {
            placeslibres[k] = emploisdisp[ ts_jr_dist[debut + k] ]; 
          }
          
          if (mode == "subjectif_d") {
            for (auto k = 0; k < k_valid; ++k) {
              attractivites[k] = ts_emplois[ ts_jr_dist[debut + k] ] * odsub[k]; 
            }
          } else if (mode== "subjectif_c") {
            for (auto k = 0; k < k_valid; ++k) {
              attractivites[k] = placeslibres[k] * odsub[k]; 
            }
          } else if (mode == "discret") {
            for (auto k = 0; k < k_valid; ++k) {
              attractivites[k] = ts_emplois[ ts_jr_dist[debut + k] ]; 
            }
          } else {
            for (auto k = 0; k < k_valid; ++k) {
              attractivites[k] = placeslibres[k]; // Cas continu de l'emploi homogène. attractivité liée à la proportion d'emplois disp.
            }
          }
          
          // Nombre d'actifs partant par freq_actif paquets.
          double actifspartant = ts_actifs[i] / freq_actifs[i];
          std::vector<double> odds(ts_xr_odds.begin() + debut, ts_xr_odds.begin() + fin);
          repartition = repartir_actifs(placeslibres, attractivites, odds, ts_fuite[i], actifspartant, SEUIL_NEWTON);
          
          // Inscription des résultats locaux dans la perspective globale.
          for(auto k = 0; k < k_valid ; ++k) {
            emploisdisp[ ts_jr_dist[debut + k] ] -= repartition[k];
            liaisons[Iboot][ debut + k ] += repartition[k];
          } // k_valid
        } // ishuf
      } // iboot
    } // max
  } // omp
} // omp nthreads

NumericVector resultat(Nx);
resultat.fill(0.0);

for (auto i = 0; i < Nx; ++i) {
  for(auto Iboot = 0; Iboot < ntr; ++Iboot) {
    resultat(i) += liaisons[Iboot][i] / Nboot;
  }
}

return List::create(_("flux") = resultat);
 }