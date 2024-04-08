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

using namespace Rcpp;

// Remarque pour les fonctions pénalités: multiplier la fonction par un facteur arbitraire ne change pas le résultat de meaps.
//' Fonction de pénalité "marche" : vaut 1 sur un rayon fixé, et decru au-delà.
 inline double marche(double x, const double rayon, const double decru) {
   if (x < rayon) return 1;
   return decru;
 }
 
 // Fonction de type logistique x -> 1 + amplitude * { exp(-(x-rayon)) - 1} / { exp(-(x-rayon)) + 1 }
 inline double logistique(double x, const double rayon, const double amplitude) {
   double ex = exp(rayon-x);
   return 1 + amplitude * ( ex - 1) / (ex + 1);
 }
 
 //' La fonction meaps en mode continu sur plusieurs shufs avec en entrée une row sparse matrix destructurée selon ses éléments.
 //' @param j_dist Le vecteur des indices des colonnes non vides.
 //' @param p_dist Le vecteur du nombres de valeurs non nulles sur chacune des lignes.
 //' @param x_dist Le vecteur des valeurs dans l'ordre de j_dist.
 //' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes). 
 //' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
 //' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude. 
 //' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i 
 //' en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs. 
 //' @param attraction Choix de la fonction de pénalité appliquée à l'accessibilité.
 //' @param alpha Premier paramètre de la fonction de pénalité, si nécessaire.
 //' @param beta Second paramètre.
 //' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto. 
 //' @param progress Ajoute une barre de progression. Default : true. 
 //' @param normalisation Calage des emplois disponibles sur le nombre d'actifs travaillant sur la zone. Default : false.
 //' @param fuite_min Seuil minimal pour la fuite d'un actif. Doit être supérieur à 0. Défault = 1e-3.
 //'
 //' @return renvoie une matrice avec les estimations du nombre de trajets de i vers j.
 // [[Rcpp::export]]
 NumericVector meaps_continu_test(IntegerVector j_dist, IntegerVector p_dist, NumericVector x_dist, NumericVector emplois,
                                 NumericVector actifs, NumericVector f, IntegerMatrix shuf, std::string attraction = "constant",
                                 double alpha = 10, double beta = 1,
                                 int nthreads = 0, bool progress = true, bool normalisation = false, double fuite_min = 1e-3) {
   const std::size_t N = actifs.size(), Nboot = shuf.nrow(), Ns = shuf.ncol();
   
   // Quelques vérifications préalables.
   // if (p_dist.size() != (N + 1L)) {
   //   stop("Le vecteur des actifs et le nb de colonnes de la matrice des distances ne correspondent pas.");
   // }
   // if (f.size() != N) {
   //   stop("Le vecteur fuite et le nb de lignes de la matrice des distances ne correspondent pas.");
   // }
   // if (actifs.size() != max(j_dist)) {
   //   stop("Le vecteur actifs et le nb de lignes de la matrice des distances ne correspondent pas (ou la dernière ligne est vide).");
   // }
   // 
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
   std::vector<int> shufcpp = as< std::vector<int> >(shuf);
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

   // Le vecteur shuf peut être plus long que le nombre de lignes de rkdist
   // s'il fait repasser plusieurs fois la même ligne d'actifs. Dans ce cas, on
   // compte la fréquence de passage de chaque ligne et l'on divise le poids de
   // la ligne par cette fréquence.
   std::vector<int> freq_actifs(N, 0L);
   for (auto i : ishuf[0]) {
     freq_actifs[i]++;
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
         // Increment progress_bar.
         p.increment();

         // Construction de l'accessibilité dite pénalisée.
         std::size_t debut = p_dist(from), fin = p_dist(from + 1L);
         std::size_t k_valid = fin - debut;

         std::vector<double> pen_attract(k_valid);  // emplois restant fois pénalité sur le site j.
         double emplois_libres = 0;
         double attract;
         for (std::size_t k = 0; k < k_valid; ++k) {
           if (attraction == "constant") {
             attract = 1;
           } else if (attraction == "marche") {
             attract = marche(xr_dist[debut + k], alpha, beta);
           } else if (attraction == "logistique") {
             attract = logistique(xr_dist[debut + k], alpha, beta);
           }
           pen_attract[k] = emp[ jr_dist[ debut + k] ] * attract;
           emplois_libres += emp[ jr_dist[ debut + k] ];
         }

         // Tests de manque de places libres.
         if (emplois_libres == 0) continue;

         double actifspartant = actifscpp[from] / freq_actifs[from];
         if (emplois_libres <= actifspartant) {
           for (std::size_t k = 0; k < k_valid; ++k) {
             liaisons[debut + k] += emp[ jr_dist[debut + k] ];
             emp[ jr_dist[debut + k] ] = 0;
           }
           continue;
         }

         std::vector<double> pen_access(k_valid);
         for (std::size_t k = 0; k < k_valid;) {
           double tot = 0;
           auto pos = k + 1L;
           while (xr_dist[debut + k] == xr_dist[debut + pos] && pos < k_valid) ++pos;
           for (std::size_t ego = k; ego < pos; ++ego) {
             tot += pen_attract[ego];
           }
           for (std::size_t ego = k; ego < pos; ++ego) {
             pen_access[ego] = tot;
           }
           k = pos;
         }

         // Calcul de l'absorption effective.
         double emploi_pen_disp = *std::max_element(pen_access.begin(), pen_access.end());
         if (emploi_pen_disp == 0) continue;  // Au cas où il n'y a pas de places libres.
         double absorption = -log(fcpp[from]) / emploi_pen_disp;

         // Calcul des actifs absorbés par sites (sans traitement des sites à
         // distances égales).
         std::vector<double> jobtakers(k_valid);
         jobtakers[0] = actifspartant * (1 - exp(-absorption * pen_access[0]));
         for (std::size_t k = 1; k < k_valid; ++k) {
           jobtakers[k] = actifspartant * exp(-absorption * pen_access[k]) - jobtakers[k - 1];
         }

         // Répartition des jobtakers dans tous les cas.
         std::vector<double> repartition(k_valid);
         for (std::size_t k = 0; k < k_valid;) {
           double tot_attraction = 0;
           double tot_jobtakers = 0;
           auto pos = k + 1L;
           while (xr_dist[debut + k] == xr_dist[debut + pos] && pos < k_valid) ++pos;
           for (std::size_t ego = 0; ego < pos; ++ego) {
             tot_attraction = +pen_attract[ego];
             tot_jobtakers = +jobtakers[ego];
           }
           if (tot_attraction > 0) {
             for (std::size_t ego = k; ego < pos; ++ego) {
               repartition[ego] = pen_attract[ego] / tot_attraction * tot_jobtakers;
             }
           }
           k = pos;
         }

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

   return actifs;
 }
 
 
 
 