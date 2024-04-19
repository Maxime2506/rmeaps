#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <algorithm>
#include <iterator>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "repartir_continu.h"

using namespace Rcpp;

// Remarque pour les fonctions pénalités: multiplier la fonction par un facteur arbitraire ne change pas le résultat de meaps.
//' Fonction de pénalité "marche" : vaut 1 sur un rayon fixé, et decru au-delà.
 //' @param x distance.
 //' @param rayon distance de la marche.
 //' @param plancher point bas après la marche.
 //' 
 //' @return un facteur d'atraction
 inline double marche(double x, const double rayon, const double plafond) {
   if (x > rayon) return 1;
   return plafond;
 }
 
 //' Fonction de pénalité "marche_liss" : vaut 1 sur un rayon fixé, et décroit au-delà,
 //' en lissant le passge pour les distances fractionaires
 //' @param x distance.
 //' @param rayon distance de la marche.
 //' @param plancher point bas après la marche.
 //' 
 //' @return un facteur d'atraction
 inline double marche_liss(double x, const double rayon, const double plafond) {
   if (x > ceil(rayon)) return 1;
   if (x <= floor(rayon)) return plafond;
   return plafond + (1 - plafond) *  (ceil(rayon) - rayon);
 }
 
 //' Fonction de pénalité "decay" : odd = 1/d^delta + plancher,
 //' @param x distance.
 //' @param delta, exposant de la distance.
 //' @param plancher valeur pour les distances infinies.
 //' 
 //' @return un facteur d'atraction
 inline double decay(double x, const double delta, const double plancher) {
   return exp(-delta * log(x)) + plancher;
 }
 
 // Fonction de type logistique x -> 1 + amplitude * { exp(-(x-rayon)) - 1} / { exp(-(x-rayon)) + 1 }
 //' @param x distance.
 //' @param rayon distance de la bascule.
 //' @param amplitude raideur de la bascule.
 //' @param plancher point bas après la marche.
 //' 
 //' @return un facteur d'atraction
 inline double logistique(double x, const double rayon, const double amplitude, const double plancher) {
   double ex = exp( (rayon-x)/amplitude );
   return plancher + ex / (ex + 1);
 }
 
 //' La fonction meaps en mode continu sur plusieurs shufs avec en entrée une Row Sparse Matrix destructurée selon ses éléments.
 //' @param j_dist Le vecteur des indices des colonnes non vides.
 //' @param p_dist Le vecteur du nombres de valeurs non nulles sur chacune des lignes.
 //' @param x_dist Le vecteur des valeurs dans l'ordre de j_dist.
 //' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes). 
 //' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
 //' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude. 
 //' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i 
 //' en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs. 
 //' @param attraction Choix de la fonction d'attraction des différents sites, appliquée à l'accessibilité. 
 //' Par défaut, "constant" où aucun site n'a plus d'attrait qu'un autre. 
 //' "marche" où l'attrait vaut 1 jusqu'à une certaine distance (param 1) puis moins (param 2). f(x) = 1 si x < p1, = p2 si x > p1.
 //' "logistique" où l'attrait décroît selon une fonction logistique avec une distance de bascule (param 1), une vitesse de bascule (param 2) 
 //' et un seuil (param p). Si h(x) = exp( (x-p1)/p2), f(x) = p3 + h(x) / (1 + h(x)).
 //' "odds" où chaque flux (from, to) se voit attribuer un odds. Dans ce cas, on entre un Row Sparse Matrix des log(odds) selon ses éléments.
 //' @param param est un vecteur avec dans l'ordre les valeurs des paramètres.
 //' @param j_odds, p_odds et x_odds sont les vecteurs de la Row Sparse Matrix lorsque attraction = "odds".
 //' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto. 
 //' @param progress Ajoute une barre de progression. Default : true. 
 //' @param normalisation Calage des emplois disponibles sur le nombre d'actifs travaillant sur la zone. Default : false.
 //' @param fuite_min Seuil minimal pour la fuite d'un actif. Doit être supérieur à 0. Défault = 1e-3.
 //'
 //' @return renvoie un vecteur des estimations des flux de i vers j.
 // [[Rcpp::export(.meaps_continu)]]
 NumericVector meaps_continu_cpp(IntegerVector j_dist, IntegerVector p_dist, NumericVector x_dist, NumericVector emplois,
                                 NumericVector actifs, NumericVector f, IntegerMatrix shuf, 
                                 NumericVector param,
                                 NumericVector j_odds,
                                 NumericVector p_odds ,
                                 NumericVector x_odds,
                                 std::string attraction = "constant",
                                 int nthreads = 0, bool progress = true, bool normalisation = false, double fuite_min = 1e-3) {
   const std::size_t N = actifs.size(), Nboot = shuf.nrow(), Ns = shuf.ncol();
   
#ifdef _OPENMP
   int ntr = nthreads;
   if (ntr == 0) {
     ntr = omp_get_max_threads();
   }
   if (ntr > omp_get_max_threads()) {
     ntr = omp_get_max_threads();
   }
   if (progress == TRUE) REprintf("Nombre de threads = %i\n", ntr);
#endif
   
   // Les rangs dans shuf commencent à 1.
   shuf = shuf - 1L;
   
   // Conversion en vecteur de vecteurs C++.
   std::vector< std::size_t > shufcpp = as< std::vector<std::size_t> >(shuf);
   std::vector< std::vector<int> > ishuf(Nboot, std::vector<int>(Ns));
   for (std::size_t i = 0; i < Nboot; ++i) {
     for (std::size_t j = 0; j < Ns; ++j) {
       ishuf[i][j] = shufcpp[i + j * Nboot];
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
   
   // REMARQUE : les conversions de Numeric et IntegerVector en équivalent std::vector sont ici nécessaires car
   // les vecteurs initiaux de sont pas thread safe dans openmp.
   // Conversion en vecteurs C++.
   std::vector<double> emploisinitial = as<std::vector<double>>(emplois);
   std::vector<double> fcpp = as<std::vector<double>>(f);
   std::vector<double> actifscpp = as<std::vector<double>>(actifs);
   
   // Passage à une Ranked Row Sparse Matrix.
   int Ndata = x_dist.size();
   std::vector<double> xr_dist(Ndata);
   std::vector<int> jr_dist(Ndata);
   
   for (std::size_t from = 0; from < N; ++from) {
     std::multimap<double, int> arrangement;
     for (std::size_t k = p_dist(from); k < p_dist(from + 1L); ++k) {
       arrangement.insert(std::make_pair(x_dist(k), j_dist(k)));
     }
     int index = 0;
     for (auto it = arrangement.begin(); it != arrangement.end(); ++it) {
       xr_dist[p_dist(from) + index] = (it->first);
       jr_dist[p_dist(from) + index] = (it->second);
       ++index;
     }
   }
   
   // paramètres pour la fonction d'attraction choisi
   std::vector<double> parametres = as< std::vector<double> >(param);
   
   std::vector<double> xr_lor;
   std::vector<int> jr_lor;
   // Le vecteur shuf peut être plus long que le nombre de lignes de rkdist
   // s'il fait repasser plusieurs fois la même ligne d'actifs. Dans ce cas, on
   // compte la fréquence de passage de chaque ligne et l'on divise le poids de
   // la ligne par cette fréquence.
   std::vector<int> freq_actifs(N, 0L);
   for (auto i : ishuf[0]) {
     freq_actifs[i]++;
   }
   
   if (attraction == "odds") {
     // Il s'agit ici de garder la propriété sparse des odds en se replaçant dans l'ordre des rangs des distances.
     for (std::size_t from = 0; from < N; ++from) {
       std::size_t odds_debut = p_odds(from), odds_fin = p_odds(from + 1L);
       for (std::size_t k = p_dist(from); k < p_dist(from + 1L); ++k) {
         auto pos = std::lower_bound(j_odds.begin() + odds_debut, j_odds.begin() + odds_fin, jr_dist[k]);
         if (pos != j_odds.begin() + odds_fin) {
           auto ind = std::distance(j_odds.begin(), pos);
           xr_lor.push_back(x_odds[ind]);
           jr_lor.push_back(j_odds[ind]);
         }
       }
     }
   }
   
   // Initialisation du résultat.
   std::vector<double> liaisons(x_dist.size(), 0.0);
   
   // Lancement du bootstrap.
   Progress p(Nboot * Ns, progress);
   
#ifdef _OPENMP
#pragma omp declare reduction(vsum : std::vector<double> : std::transform(                     \
   omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>()))      \
     initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp parallel for num_threads(ntr) \
     shared(Nboot, N, ishuf, emploisinitial, fcpp, actifscpp, xr_dist, jr_dist, p_dist) reduction(vsum : liaisons)
#endif
     
     for (int iboot = 0; iboot < Nboot; ++iboot) {
       // Initialisation de de l'ordre de départ des actifs et de l'emploi
       // disponible au début.
       std::vector<int> theshuf = ishuf[iboot];  // deep copy pour un boot.
       std::vector<double> emp(emploisinitial);  // deep copy.
       
       for (auto from : theshuf) {
         // check interrupt & progress
         if (! Progress::check_abort() )
           p.increment();
         // Construction de l'accessibilité dite pénalisée.
         std::size_t debut = p_dist(from), fin = p_dist(from + 1L);
         std::size_t k_valid = fin - debut;
         
         std::vector<double> facteur_attraction(k_valid, 1.0), emplois_libres(k_valid), repartition(k_valid);
         std::size_t odds_index = 0;
         
         if (attraction == "marche") {
           for (std::size_t k = 0; k < k_valid; ++k) {
             facteur_attraction[k]  = marche(xr_dist[debut + k], parametres[0], parametres[1]);
           }}
         
         if (attraction == "marche_liss") {
           for (std::size_t k = 0; k < k_valid; ++k) {
             facteur_attraction[k]  = marche_liss(xr_dist[debut + k], parametres[0], parametres[1]);
           }}
         
         if (attraction == "decay") {
           for (std::size_t k = 0; k < k_valid; ++k) {
             facteur_attraction[k]  = decay( xr_dist[debut + k], parametres[0], parametres[1]);
           }}
         
         if (attraction == "logistique") {
           for (std::size_t k = 0; k < k_valid; ++k) {
             facteur_attraction[k] = logistique(xr_dist[debut + k], parametres[0], parametres[1], parametres[2]);
           }}
         
         if (attraction == "odds") {
           std::size_t odds_index = 0;
           for (std::size_t k = 0; k < k_valid; ++k) {
             if (jr_lor[p_odds(from) + odds_index] == jr_dist[debut + k]) {
               facteur_attraction[k] = exp( xr_lor[p_odds(from) + odds_index] );
               odds_index++;
             } 
           }
         }
         
         for (std::size_t k = 0; k < k_valid; ++k) {
           emplois_libres[k] = emp[ jr_dist[ debut + k] ];
         } 
         
         double actifspartant = actifscpp[from] / freq_actifs[from];
         
         std::vector<double> dist(xr_dist.begin() + debut, xr_dist.begin() + fin);
         repartition = repartir_continu(actifspartant, fcpp[from], facteur_attraction, dist, emplois_libres); 
         
         // Impact sur l'emploi disponible total et sommation sur les emplois pris.
         for (std::size_t k = 0; k < k_valid; ++k) {
           emp[ jr_dist[debut + k] ] -= repartition[k];
           liaisons[debut + k] += repartition[k];
         }
       }
     }
     
     // Il faut ressortir le résultat au format classic de Sparse Matrix en
     // réordonnant les valeurs par l'ordre des colonnes.
     NumericVector resultat(liaisons.size());
   for (std::size_t from = 0; from < N; ++from) {
     int debut = p_dist[from], 
                       fin = p_dist[from + 1L];
     for (std::size_t k = debut; k < fin; ++k) {
       auto pos = std::lower_bound(j_dist.begin() + debut, j_dist.begin() + fin, jr_dist[k]);
       auto ind = std::distance(j_dist.begin(), pos);
       resultat[ind] = liaisons[k] / Nboot;
     }
   }
   
   return resultat;
 }
 
 
 
 
 
 