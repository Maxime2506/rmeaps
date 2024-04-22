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
 // [[Rcpp::export(.another_meaps)]]
 Rcpp::DataFrame another_meaps_cpp(IntegerVector jr_dist, 
                               IntegerVector p_dist, 
                               NumericVector xr_dist, 
                               NumericVector emplois,
                               NumericVector actifs, 
                               NumericVector f, 
                               NumericVector param,
                               IntegerVector jr_odds,
                               IntegerVector p_odds ,
                               NumericVector xr_odds,
                               std::string attraction = "constant",
                               int nthreads = 0, bool progress = true, bool normalisation = false, double fuite_min = 1e-3) {
   const int _LIMITE_LOOP = 1000;
   const double _LIMITE_PRECISION = 1e-9;
   
   const std::size_t N = actifs.size(), K = emplois.size();
   
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
   
   // Choix d'une limite basse pour la fuite.
   f = ifelse(f > fuite_min, f, fuite_min);
   
   // Calage de l'emploi sur les actifs.
   if (normalisation) {
     emplois = emplois * sum(actifs * (1 - f)) / sum(emplois);
   }
   
   std::vector<double> emplois_libres = as<std::vector<double>>(emplois);
   std::vector<double> fcpp = as<std::vector<double>>(f);
   std::vector<double> actifs_libres = as<std::vector<double>>(actifs);
  
   std::vector<double> parametres = as< std::vector<double> >(param);
   
   std::vector<int> _jr_dist = as< std::vector<int> >(jr_dist);
   std::vector<int> _p_dist = as< std::vector<int> >(p_dist);
   std::vector<double> _xr_dist = as< std::vector<double> >(xr_dist);
  
   std::vector<int> _jr_odds = as< std::vector<int> >(jr_odds);
   std::vector<int> _p_odds = as< std::vector<int> >(p_odds);
   std::vector<double> _xr_odds = as< std::vector<double> >(xr_odds);
   
   // Initialisation du résultat.
   std::vector< std::vector<double> > liaisons(N, std::vector<double> (K));
   
   double tot_actifs_libres = 0, old_tot;
   for (int i = 0; i < N; ++i) {
     tot_actifs_libres += actifs_libres[i];
   }
   
   int nloop = 0;
   
   do {
   
   // Lancement du bootstrap.
   if (progress == TRUE) REprintf("boucle n° %i\n", nloop + 1L);
#ifdef _OPENMP
#pragma omp parallel for num_threads(ntr)  
#endif
       
       for (std::size_t from = 0; from < N; ++from) {
         
         // Curseurs de la ligne from en cours.
         std::size_t debut = _p_dist[from], fin = _p_dist[from + 1L];
         std::size_t k_valid = fin - debut;
         
         // Calcul de la fonction d'attraction retenue et mise en rang des emplois cibles.
         std::vector<double> facteur_attraction(k_valid), cibles(k_valid), repartition(k_valid);
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
           cibles[k] = emplois_libres[ _jr_dist[ debut + k] ];
         }
         
         
         // Calcul de la distribution sur la ligne.
         std::vector<double> dist(_xr_dist.begin() + debut, _xr_dist.begin() + fin);
         repartition = one_distrib_continu(actifs_libres[from], fcpp[from], facteur_attraction, dist, cibles);
         
         // Réagencement des résultats.
         for (std::size_t k = 0; k < k_valid; ++k) {
           liaisons[from][_jr_dist[debut + k]] = repartition[k];
         }
       } 
     
     // Après la distribution simultanée, il faut récupérer les excédents et reconstruire des vecteurs
     // d'actifs encore libres et d'emplois encore libres.
     for (std::size_t i = 0; i < N; ++i) {
       actifs_libres[i] = 0;
     }
     double tx_depassement;
     // Traitement des dépassements en colonnes.
     for (std::size_t j = 0; j < K; ++j) {
       emplois_libres[j] = emplois[j];
       for (std::size_t i = 0; i < N; ++i) {
         emplois_libres[j] -= liaisons[i][j];
       }
       if (emplois_libres[j] < 0) {
         tx_depassement = emplois[j] / (emplois[j] - emplois_libres[j]);
         // Renvoie à domicile des actifs excédentaires à proportion de leur contribution à l'excédent.
         // ATTENTION : le renvoi à domicile est à proportion du total des actifs en place, et non pas 
         // uniquement à proportion de la dernière vague d'arrivants.
         for (std::size_t i = 0; i < N; ++i) {
           actifs_libres[i] += liaisons[i][j] * (1 - tx_depassement);
           liaisons[i][j] = liaisons[i][j] * tx_depassement; 
           }
         emplois_libres[j] = 0;
       }
      }
       
     
     old_tot = tot_actifs_libres;
     tot_actifs_libres = 0;
     for (std::size_t i = 0; i < N; ++i) {
       actifs_libres[i] /= (1 - fcpp[i]); // Ne pas oublier de renvoyer aussi à domicile les fuyards des actifs renvoyés.
       tot_actifs_libres += actifs_libres[i];
     }
     nloop++;
     
   } while (tot_actifs_libres > 0 && std::abs(tot_actifs_libres - old_tot) > 0.1 && nloop < _LIMITE_LOOP);
     
     
   std::vector<int> i_res, j_res;
   std::vector<double> x_res;
   for (std::size_t i = 0; i < N; ++i) {
     for (std::size_t j = 0; j < K; ++j) {
       if (liaisons[i][j] > _LIMITE_PRECISION) {
         i_res.push_back(i + 1L);
         j_res.push_back(j + 1L);
         x_res.push_back(liaisons[i][j]);
       }
     }
   }
   
   return Rcpp::DataFrame::create(_("i") = i_res, _("j") = j_res, _("flux") = x_res);
 }
 
 
 
 