#ifdef _OPENMP
# include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "repartir_actifs.h"
#include "calculer_rang.h"

using namespace Rcpp;

 //' La fonction meaps sur plusieurs shufs
 //' @param dist La matrice (creuses en colonnes) des distances dans lequel les colonnes j sont passées en revue pour chacune des lignes i.
 //' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes).
 //' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
 //' @param modds La matrice des odds modifiant la chance d'absorption de chacun des sites j pour des résidents en i.
 //' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude.
 //' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs.
 //' @param mode Choix du rôle de l'emploi disponible au cours du processus. Default : continu. Autre choix : discret, subjectif_c ou _d...
 //' @param odds_subjectifs Attractivité des sites proches des actifs, pour le mode defini. default : null.
 //' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto.
 //' @param progress Ajoute une barre de progression. Default : true.
 //' @param normalisation Calage des emplois disponibles sur le nombre d'actifs travaillant sur la zone. Default : false.
 //' @param fuite_min Seuil minimal pour la fuite d'un actif. Doit être supérieur à 0. Défault = 1e-3.
 //' 
 //' @return renvoie une matrice avec les estimations du nombre de trajets de i vers j.
 // [[Rcpp::export]]
 NumericVector meaps_dgr2(IntegerVector j_dist,
                         IntegerVector p_dist,
                         NumericVector x_dist,
                         const int N,
                         const int K,
                         NumericVector emplois,
                         NumericVector actifs,
                         IntegerVector j_modds,
                         IntegerVector p_modds,
                         NumericVector x_modds,
                         NumericVector f,
                         IntegerMatrix shuf,
                         std::string mode = "continu",
                         Nullable<NumericVector> oddssubjectifs = R_NilValue,
                         int nthreads = 0,
                         bool progress = true,
                         bool normalisation = false,
                         double fuite_min = 1e-3,
                         double seuil_newton = 1e-6) {
   
   const int Nboot = shuf.nrow(),
     Ns = shuf.ncol();
   
   // Quelques vérifications préalables.
   if (emplois.size() != K) {
     stop("Le vecteur emplois et le nb de colonnes de la matrice rkdist ne correspondent pas.");
   }
   if (f.size() != N) {
     stop("Le vecteur fuite et le nb de lignes de la matrice rkdist ne correspondent pas.");
   }
   if (actifs.size() != N) {
     stop("Le vecteur actifs et le nb de lignes de la matrice rkdist ne correspondent pas.");
   }
   
   if ((mode == "subjectif_c" || mode == "subjectif_d")&& oddssubjectifs.isNull()) {
     stop("En mode subjectif, le vecteur oddssubjectifs doit être défini");
   }
   
   std::vector<double> odsub;
   if (mode == "subjectif_c"|| mode == "subjectif_d") {
     odsub = as< std::vector<double> >(oddssubjectifs);
   }
   
#ifdef _OPENMP
   int ntr = nthreads;
   if (ntr == 0) {ntr = omp_get_max_threads();}
   if (ntr > omp_get_max_threads()) { ntr = omp_get_max_threads(); }
   if (progress==TRUE) REprintf("Nombre de threads = %i\n", ntr);
#endif
   
   // Les rangs dans shuf commencent à 1.
   shuf = shuf - 1L;
   
   // Conversion en vecteur de vecteurs C++.
   std::vector<int> shufcpp = as< std::vector<int> >(shuf);
   std::vector< std::vector<int> > ishuf(Nboot, std::vector<int>(Ns));
   for (int i = 0; i < Nboot; ++i) {
     for (int j = 0; j < Ns; ++j) {
       ishuf[i][j] = shufcpp[ i  + j * Nboot];
       if (ishuf[i][j] >= N) { 
         stop("Les rangs de shufs vont au-delà de la taille du vecteur actifs.");
       }
     }
   }
   
   // Choix d'une limite basse pour la fuite.
   f = ifelse(f > fuite_min, f, fuite_min);
   
   // Calage de l'emploi sur les actifs.
   if (normalisation) {
     emplois = emplois * sum(actifs * (1 - f)) / sum(emplois); 
   }
   
   // Conversion en vecteurs c++ 
   std::vector<double> emploisinitial = as< std::vector<double> >(emplois);
   std::vector<double> fcpp = as< std::vector<double> >(f);
   std::vector<double> actifscpp = as< std::vector<double> >(actifs);
   
   // Préparation des données @x de odds sur la structure @p @j des distances.
   std::vector<double> odds_x(x_dist.size());
   for (int from = 0; from < N; ++from) {
     for (int k = p_dist(from); k != p_dist(from + 1L); ++k) {
       auto pos = std::find(j_modds.begin() + p_modds(from), j_modds.begin() + p_modds(from + 1L), j_dist[k]);
       if (*pos == j_modds[p_modds(from + 1L)]) {
         odds_x.push_back(1);
       } else {
         double val = exp( x_modds[*pos]);
         odds_x.push_back(val);
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
   
   // Initialisation du résultat.
   // On conserve l'ordre subjectif et on ne remplit que les non zeros (différent de la version précédente)
   // std::vector< std::vector<double> > liaisons(N, std::vector<double>(K));
   // for (int from = 0; from < N; ++i) {
   //   liaisons[from].resize(dist.outerIndexPtr(from + 1L) - dist.outerIndexPtr(from))
   // }
   
   std::vector<double> liaisons(x_dist.size(), 0.0);
   
   // Lancement du bootstrap.
   Progress p(Nboot * Ns, progress);
#ifdef _OPENMP
#pragma omp declare reduction(vsum : std::vector<double> :               \
   std::transform(omp_out.begin(), omp_out.end(),                        \
                  omp_in.begin(), omp_out.begin(), std::plus<double>())) \
     initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp parallel for num_threads(ntr)                                                          \
     shared(Nboot, N, ishuf, emploisinitial, fcpp, actifscpp, x_dist, p_dist, odds_x)                                \
     reduction (vsum : liaisons)
#endif
     for (int iboot = 0; iboot < Nboot; ++iboot) {
       
       std::vector<int> theshuf = ishuf[iboot]; // deep copy pour un boot.
       
       // Initialisation de l'emploi disponible.
       std::vector<double> emp(emploisinitial); // deep copy.
       
       for (auto from: theshuf) {
         
         // Increment progress_bar.
         p.increment();
         
         // Construction des vecteurs dans l'ordre des rangs de distances
         int k_valid = p_dist(from + 1L) - p_dist(from);
         std::vector<int> arr(k_valid);
         std::vector<double> les_distances(x_dist.begin() + p_dist(from), x_dist.begin() + p_dist(from + 1L));
         
         arr = calculer_rang(les_distances);
         
         std::vector<int> index(k_valid);
         std::iota(begin(index), end(index), 0);
         std::sort(begin(index), end(index), 
                   [&arr](std::size_t i, std::size_t j) { return arr[i] < arr[j]; });
         
         std::vector<double> odds(k_valid);
         for (int k = 0; k < k_valid; ++k) {
           odds[k] = odds_x[ p_dist(from) + index[k] ];
         }
         
         std::vector<double> placeslibres (k_valid), attractivites(k_valid), repartition(k_valid);
         for (int k = 0; k < k_valid; ++k) {
           placeslibres[k] = emp[ arr[k] ]; 
         }
         
         if (mode == "subjectif_d") {
           for (int k = 0; k < k_valid; ++k) {
             attractivites[k] = emploisinitial[ arr[k] ] * odsub[k]; 
           }
         } else if (mode== "subjectif_c") {
           for (int k = 0; k < k_valid; ++k) {
             attractivites[k] = placeslibres[k] * odsub[k]; 
           }
         } else if (mode == "discret") {
           for (int k = 0; k < k_valid; ++k) {
             attractivites[k] = emploisinitial[ arr[k] ]; 
           }
         } else {
           for (int k = 0; k < k_valid; ++k) {
             attractivites[k] = placeslibres[k]; // Cas continu de l'emploi homogène. attractivité liée à la proportion d'emplois disp.
           }
         }
         
         // Nombre d'actifs partant par freq_actif paquets.
         double actifspartant = actifscpp[from] / freq_actifs[from];
         
         repartition = repartir_actifs(placeslibres, attractivites, odds, fcpp[from], actifspartant, seuil_newton);
         
         // Inscription des résultats locaux dans la perspective globale.
         std::vector<int> indexes = arr;
         std::sort(indexes.begin(), indexes.end());
         
         for(int k = 0; k < k_valid ; ++k) {
           emp[ arr[k] ] -= repartition[k];
           
           int indice = *std::find(indexes.begin(), indexes.end(), arr[k]); //Nécessaire car arr change à chaque boot
           liaisons[ p_dist[from] + indice ] += repartition[k]; 
         }
       }
     }
 
   NumericVector resultat = wrap(liaisons);
   
   // division par le nb de boot.
   resultat = resultat / Nboot ;
   
   return resultat;
 }
 
 
 
 