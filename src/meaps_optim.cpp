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
 
 //' La fonction meaps_continu qui ne renvoit que le KL de l'estimation en référence à une distribution connue. 
 //' @param jr_dist Le vecteur des indices des colonnes non vides.
 //' @param p_dist Le vecteur du nombres de valeurs non nulles sur chacune des lignes.
 //' @param xr_dist Le vecteur des valeurs dans l'ordre de jr_dist.
 //' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes). 
 //' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
 //' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude. 
 //' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. Il est possible de segmenter les départs d'une ligne i 
 //' en répétant cette ligne à plusieurs endroits du shuf et en répartissant les poids au sein du vecteurs actifs. 
 //' @param row_group Le vecteur de regroupement des départs (par ex. code commune en integer).
 //' @param col_group Le vecteur de regroupement des arrivées (par ex. code commune en integer).
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
 //' @return renvoie une matrice avec les estimations des flux regroupés.
 // [[Rcpp::export(.meaps_optim)]]
 NumericMatrix meaps_optim_cpp(IntegerVector jr_dist, 
                               IntegerVector p_dist, 
                               NumericVector xr_dist, 
                               NumericVector emplois,
                               NumericVector actifs, 
                               NumericVector f, 
                               IntegerMatrix shuf, 
                               IntegerVector row_group,
                               IntegerVector col_group,
                               NumericVector param,
                               IntegerVector jr_odds,
                               IntegerVector p_odds ,
                               NumericVector xr_odds,
                               std::string attraction = "constant",
                               int nthreads = 0, bool progress = true, bool normalisation = false, double fuite_min = 1e-3) {
  const std::size_t N = actifs.size(), Nboot = shuf.nrow(), Ns = shuf.ncol();

auto Nref = *std::max_element(row_group.begin(), row_group.end());
Nref = Nref + 1L;
auto Kref = *std::max_element(col_group.begin(), col_group.end());
Kref = Kref + 1L;

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

   // REMARQUE : les conversions de Numeric et IntegerVector en équivalent std::vector sont ici nécessaires car
   // les vecteurs initiaux de sont pas thread safe dans openmp.
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

   std::vector<double> emploisinitial = as<std::vector<double>>(emplois);
   std::vector<double> fcpp = as<std::vector<double>>(f);
   std::vector<double> actifscpp = as<std::vector<double>>(actifs);

   std::vector<double> parametres = as< std::vector<double> >(param);

   std::vector<int> _jr_dist = as< std::vector<int> >(jr_dist);
   std::vector<int> _p_dist = as< std::vector<int> >(p_dist);
   std::vector<double> _xr_dist = as< std::vector<double> >(xr_dist);

   std::vector<int> _jr_odds = as< std::vector<int> >(jr_odds);
   std::vector<int> _p_odds = as< std::vector<int> >(p_odds);
   std::vector<double> _xr_odds = as< std::vector<double> >(xr_odds);


   // Choix d'une limite basse pour la fuite.
   f = ifelse(f > fuite_min, f, fuite_min);

   // Calage de l'emploi sur les actifs.
   if (normalisation) {
     emplois = emplois * sum(actifs * (1 - f)) / sum(emplois);
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
   // Un vecteur représentant la matrice des flux groupés.
   std::vector<double> liaisons(Nref * Kref, 0.0);

   // Lancement du bootstrap.
   Progress p(Nboot * Ns, progress);
#ifdef _OPENMP
#pragma omp declare reduction(vsum : std::vector<double> : std::transform(                     \
   omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>()))      \
     initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp parallel for num_threads(ntr) \
     shared(Nboot, N, ishuf, emploisinitial, fcpp, actifscpp, _xr_dist, _jr_dist, _p_dist, _xr_odds, _p_odds, _jr_odds) reduction(vsum : liaisons)
#endif

     for (int iboot = 0; iboot < Nboot; ++iboot) {
       // Initialisation de l'ordre de départ des actifs et de l'emploi
       // disponible au début.
       std::vector<int> theshuf = ishuf[iboot];  // deep copy pour un boot.
       std::vector<double> emp(emploisinitial);  // deep copy.

       for (auto from : theshuf) {
         // Increment progress_bar.
         p.increment();

         // Construction de l'accessibilité dite pénalisée.
         std::size_t debut = _p_dist[from], fin = _p_dist[from + 1L];
         std::size_t k_valid = fin - debut;

         std::vector<double> facteur_attraction(k_valid), emplois_libres(k_valid), repartition(k_valid);
         std::size_t odds_index = 0;
         for (std::size_t k = 0; k < k_valid; ++k) {
           if (attraction == "constant") {
             facteur_attraction[k]  = 1;
           } else if (attraction == "marche") {
             facteur_attraction[k]  = marche(_xr_dist[debut + k], parametres[0], parametres[1]);
           } else if (attraction == "logistique") {
             facteur_attraction[k] = logistique(_xr_dist[debut + k], parametres[0], parametres[1], parametres[2]);
           } else if (attraction == "odds") {
             if (_jr_odds[_p_odds[from] + odds_index] == _jr_dist[debut + k]) {
               facteur_attraction[k] = exp( _xr_odds[_p_odds[from] + odds_index] );
               odds_index++;
             } else {
               facteur_attraction[k] = 1;
             }
           }
           emplois_libres[k] = emp[ _jr_dist[ debut + k] ];
         }

         double actifspartant = actifscpp[from] / freq_actifs[from];

         std::vector<double> dist(_xr_dist.begin() + debut, _xr_dist.begin() + fin);
         repartition = repartir_continu(actifspartant, fcpp[from], facteur_attraction, dist, emplois_libres);

         // Impact sur l'emploi disponible total et sommation sur les emplois pris.
         std::size_t curseur_ligne = row_group[from] * Nref;
         for (std::size_t k = 0; k < k_valid; ++k) {
           emp[ _jr_dist[debut + k] ] -= repartition[k];

           liaisons[ curseur_ligne + col_group[_jr_dist[debut + k]] ] += repartition[k];
         }
       }
     }

     NumericMatrix resultat(Nref, Kref);
     for (std::size_t i = 0; i < Nref; ++i) {
       for (std::size_t j = 0; j < Kref; ++j) {
         resultat(i,j) = liaisons[Nref * i + j] / Nboot;
       }
     }

     return resultat;
 }


 // Métrique pour comparer les flux agrégés estimés et des flux cibles.
 // On calcule plus simplement sur les effectifs que sur les probabilités. Ne fait que translater le résultat.
 // Si on note Ne = sum(estim) et Nc = sum(cible), on a :
 // objectif_kl = Ne.KL + Ne.log(Ne/Nc)
 // Attention : les flux estimés f_ij qui n'ont pas de flux cible correspondant (cible_ij = 0) sortent du calcul.
 double objectif_kl (NumericMatrix estim, NumericMatrix cible, double pseudozero = 1e-3) {

   if (estim.nrow() != cible.nrow()) stop("Pb sur le nombre de lignes agrégées.");
   if (estim.ncol() != cible.ncol()) stop("Pb sur le nombre de colonnes agrégées");

     double tot = 0;
     for (auto i = 0; i < cible.nrow(); ++i)
       for (auto j = 0; j < cible.ncol(); ++j) {
         if (cible(i,j) != 0 && estim(i,j) != 0) tot += estim(i,j) * (log(estim(i,j)) - log(cible(i,j)));
       }
   return tot;
 }




