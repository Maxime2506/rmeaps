#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#include <algorithm>
#include <iterator>

#include "repartir_continu.h"
#include "fcts_penal.h"

using namespace Rcpp;

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
                                   int nthreads = 0, bool verbose = true, bool normalisation = false, double fuite_min = 1e-3) {
   const int _LIMITE_LOOP = 1000;
   const double _LIMITE_PRECISION = 1e-3;
   
   const std::size_t N = actifs.size(), K = emplois.size(), Ndata = xr_dist.size();
   
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
   
   // Choix d'une limite basse pour la fuite.
   f = ifelse(f > fuite_min, f, fuite_min);
   
   // Calage de l'emploi sur les actifs.
   if (normalisation) {
     emplois = emplois * sum(actifs * (1 - f)) / sum(emplois);
   }
   
   std::vector<double> emplois_libres = as<std::vector<double>>(emplois);
   std::vector<double> _emplois = as<std::vector<double>>(emplois);
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
   
       #ifdef _OPENMP
       #pragma omp parallel for num_threads(ntr)  
       #endif
       for (std::size_t from = 0; from < N; ++from) {
         
         // Curseurs de la ligne from en cours.
         std::size_t debut = _p_dist[from], fin = _p_dist[from + 1L];
         std::size_t k_valid = fin - debut;
         
         // Calcul de la fonction d'attraction retenue et mise en rang des emplois cibles.
         std::vector<double> facteur_attraction(k_valid, 1.0), cibles(k_valid), repartition(k_valid);
         
         if (attraction == "marche") {
           for (std::size_t k = 0; k < k_valid; ++k) {
             facteur_attraction[k]  = marche(_xr_dist[debut + k], parametres[0], parametres[1]);
           }}
         
         if (attraction == "marche_liss") {
           for (std::size_t k = 0; k < k_valid; ++k) {
             facteur_attraction[k]  = marche_liss(_xr_dist[debut + k], parametres[0], parametres[1]);
           }}
         
         if (attraction == "double_marche_liss") {
           for (std::size_t k = 0; k < k_valid; ++k) {
             facteur_attraction[k]  = marche_liss(_xr_dist[debut + k], parametres[0], parametres[1], parametres[2], parametres[3]);
           }}
         
         if (attraction == "decay") {
           for (std::size_t k = 0; k < k_valid; ++k) {
             facteur_attraction[k]  = decay(_xr_dist[debut + k], parametres[0], parametres[1]);
           }}
         
         if (attraction == "logistique") {
           for (std::size_t k = 0; k < k_valid; ++k) {
             facteur_attraction[k] = logistique(_xr_dist[debut + k], parametres[0], parametres[1], parametres[2]);
           }}
         
         if (attraction == "odds") {
           std::size_t odds_index = 0;
           for (std::size_t k = 0; k < k_valid; ++k) {
             if (_jr_odds[_p_odds[from] + odds_index] == _jr_dist[debut + k]) {
               facteur_attraction[k] = exp( _xr_odds[_p_odds[from] + odds_index] );
               odds_index++;
             } 
           } 
         }
         
         for (std::size_t k = 0; k < k_valid; ++k) {
           cibles[k] = emplois_libres[ _jr_dist[ debut + k] ];
         }
         // Calcul de la distribution sur la ligne.
         std::vector<double> dist(_xr_dist.begin() + debut, _xr_dist.begin() + fin);
         repartition = one_distrib_continu(actifs_libres[from], fcpp[from], facteur_attraction, dist, cibles);
         
         // Réagencement au sein de la matrice des résultats (non sparse pour mener des calculs simples en colonnes).
         for (std::size_t k = 0; k < k_valid; ++k) {
           liaisons[from][ _jr_dist[debut + k] ] += repartition[k];
         }
       } 
     
     // Après la distribution simultanée, il faut récupérer les excédents et reconstruire des vecteurs
     // d'actifs encore libres et d'emplois encore libres.
     for (std::size_t i = 0; i < N; ++i) {
       actifs_libres[i] = 0;
     }
     double tx_depassement;
     
     // Traitement des dépassements en colonnes.
     #ifdef _OPENMP
#pragma omp declare reduction(vsum : std::vector<double> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp parallel for num_threads(ntr) reduction(vsum: actifs_libres)
     #endif
     for (std::size_t j = 0; j < K; ++j) {
       emplois_libres[j] = _emplois[j];
       for (std::size_t i = 0; i < N; ++i) {
         emplois_libres[j] -= liaisons[i][j];
       }
       if (emplois_libres[j] < 0) {
         tx_depassement = _emplois[j] / (_emplois[j] - emplois_libres[j]);

         // Renvoie à domicile des actifs excédentaires à proportion de leur contribution à l'excédent.
         // CHOIX METHODO : le renvoi à domicile est à proportion du total des actifs en place, et non pas 
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
     
   } while (tot_actifs_libres > 0 && std::abs(tot_actifs_libres - old_tot) > _LIMITE_PRECISION && nloop < _LIMITE_LOOP);
  
  
   if (nloop == _LIMITE_LOOP) REprintf("Warning : limite loop atteinte!");
   if (verbose == TRUE) {
     REprintf("Nombre de boucles effectuées = %i\n", nloop); 
     REprintf("Nombre d'actifs non occupés = %f\n", tot_actifs_libres);
     }
   
   std::vector<int> i_res(Ndata), j_res(Ndata);
   std::vector<double> x_res(Ndata);
   #ifdef _OPENMP
   #pragma omp parallel for num_threads(ntr)  
   #endif
   for (std::size_t i = 0; i < N; ++i) {
     for (std::size_t j = _p_dist[i]; j < _p_dist[i + 1L]; ++j) {
       i_res[j] = i + 1L;
       j_res[j] = _jr_dist[j] + 1L;
       x_res[j] = liaisons[i][ _jr_dist[j] ];
     }
   }
   
   return Rcpp::DataFrame::create(_("i") = i_res, _("j") = j_res, _("flux") = x_res);
 }
 
 
 
 