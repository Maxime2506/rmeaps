#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <algorithm>
#include <iterator>

#include "fcts_penal.h"
#include "meaps_core.h"

using namespace Rcpp;

 //' La fonction meaps qui distribue tous les actifs en même temps. En entrée, la matrice des distances (et si besoin des odds)
 //' doit être définie sous forme des inner et outer index d'une matrice sparse en ligne. Ceci revient à classer un data.frame
 //' avec les colonnes i, j et dist, d'abord par i, puis dist (=xr), puis j (=jr). 
 //' @param jr_dist Le vecteur des indices des colonnes non vides.
 //' @param p_dist Le vecteur du nombres de valeurs non nulles sur chacune des lignes.
 //' @param xr_dist Le vecteur des valeurs dans l'ordre de jr_dist.
 //' @param group_from Le vecteur des regroupements des lignes d'actifs.
 //' @param group_to Le vecteur des regroupements des colonnes d'emplois.
 //' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes). 
 //' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf.
 //' @param fuite Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude. 
 //' @param parametres Un vecteur avec les paramètres nécessaires selon la fonction d'attraction retenue;
 //' @param attraction Choix de la fonction d'attraction des différents sites, appliquée à l'accessibilité. 
 //' Par défaut, "constant" où aucun site n'a plus d'attrait qu'un autre. 
 //' "marche" où l'attrait vaut 1 jusqu'à une certaine distance (param 1) puis moins (param 2). f(x) = 1 si x < p1, = p2 si x > p1.
 //' "logistique" où l'attrait décroît selon une fonction logistique avec une distance de bascule (param 1), une vitesse de bascule (param 2) 
 //' et un seuil (param p). Si h(x) = exp( (x-p1)/p2), f(x) = p3 + h(x) / (1 + h(x)).
 //' "odds" où chaque flux (from, to) se voit attribuer un odds. 
 //' @param param est un vecteur avec dans l'ordre les valeurs des paramètres.
 //' @param xr_odds donne la valeur des odds en suivant la structure de jr_dist et p_dist.
 //' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto. 
 //' @param verbose Défaut = true. 
 //' @param normalisation Calage des emplois disponibles sur le nombre d'actifs travaillant sur la zone. Défaut : false.
 //' @param fuite_min Seuil minimal pour la fuite d'un actif. Doit être supérieur à 0. Défaut = 1e-3.
 //'
 //' @return renvoie les flux au format triplet.
 // [[Rcpp::export]]
 Rcpp::DataFrame meaps_all_in_optim(const IntegerVector jr_dist, 
                                    const IntegerVector p_dist, 
                                    const NumericVector xr_dist, 
                                    const IntegerVector group_from,
                                    const IntegerVector group_to,
                                    NumericVector emplois,
                                    const NumericVector actifs, 
                                    NumericVector fuite, 
                                    const NumericVector parametres,
                                    const NumericVector xr_odds,
                                    const std::string attraction = "constant",
                                    const int nthreads = 0, 
                                    const bool verbose = true, 
                                    bool normalisation = false, double fuite_min = 1e-3) {
 
   const std::size_t N = actifs.size(), K = emplois.size();
   auto Nref = 1L + *std::max_element(group_from.begin(), group_from.end());
   //Nref = Nref + 1L;
   auto Kref = 1L + *std::max_element(group_to.begin(), group_to.end());
   //Kref = Kref + 1L;

   // Choix d'une limite basse pour la fuite.
   fuite = ifelse(fuite > fuite_min, fuite, fuite_min);
   
   // Calage de l'emploi sur les actifs.
   if (normalisation) {
     emplois = emplois * sum(actifs * (1 - fuite)) / sum(emplois);
   }
   
   // Passage explicite en std::vector pour rendre les vecteurs thread safe (ts_)(nécessaire pour openmp dans la macro).
   std::vector<double> ts_emplois = as<std::vector<double>>(emplois);
   std::vector<double> ts_fuite = as<std::vector<double>>(fuite);
   const std::vector<double> ts_actifs = as<std::vector<double>>(actifs);
   
   const std::vector<double> ts_parametres = as< std::vector<double> >(parametres);
   
   const std::vector<int> ts_jr_dist = as< std::vector<int> >(jr_dist);
   const std::vector<int> ts_p_dist = as< std::vector<int> >(p_dist);
   const std::vector<double> ts_xr_dist = as< std::vector<double> >(xr_dist);
   
   const std::vector<double> ts_xr_odds = as< std::vector<double> >(xr_odds);
   
   // Initialisation du résultat.
   std::vector< std::vector<double> > liaisons(N, std::vector<double> (K));
   
   liaisons = meaps_core(ts_jr_dist, ts_p_dist, ts_xr_dist, ts_emplois, ts_actifs, ts_fuite, 
                         ts_parametres, ts_xr_odds, attraction, nthreads, verbose);

  // sortie du résultat agrégé.
  NumericMatrix agregat(Nref, Kref);
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < K; ++j) {
      agregat(group_from[i], group_from[j]) += liaisons[i][j];
    }
  }

return agregat;
 }
 
 
 
 