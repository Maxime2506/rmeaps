#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <stdio.h>
#include <algorithm>

#include "constants.h"
#include "classes.h"
#include "classes_attraction.cpp"

using namespace Rcpp;

//' La fonction meaps qui distribue tous les actifs en même temps. En entrée, la matrice des distances (et si besoin des odds)
 //' doit être définie sous forme des inner et outer index d'une matrice sparse en ligne. Ceci revient à classer un data.frame
 //' avec les colonnes i, j et dist, d'abord par i, puis dist (=xr), puis j (=jr). 
 //' @param jr_dist Le vecteur des indices des colonnes non vides.
//' @param p_dist Le vecteur du nombres de valeurs non nulles sur chacune des lignes.
 //' @param xr_dist Le vecteur des valeurs dans l'ordre de jr_dist.
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
//' @param group_from Le vecteur des regroupements des lignes d'actifs.
//' @param group_to Le vecteur des regroupements des colonnes d'emplois.
//' @param cible Le vecteur des flux groupés observés dans l'ordre des index c++ de (group_from, group_to).
//' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto. 
 //' @param verbose Défaut = true. 
//' On doit gérer normalisation et fuite_min auparavant.
 //'
//' @return renvoie les flux au format triplet.
 // [[Rcpp::export]]
List meaps(const IntegerVector jr_dist, 
                  const IntegerVector p_dist, 
                  const NumericVector xr_dist, 
                  const NumericVector emplois,
                  const NumericVector actifs, 
                  const NumericVector fuites, 
                  const NumericVector parametres,
                  const NumericVector xr_odds,
                  const std::string attraction = "constant",
                  const Nullable<IntegerVector> group_from = R_NilValue,
                  const Nullable<IntegerVector> group_to = R_NilValue,
                  const Nullable<NumericVector> cible = R_NilValue,
                  const int nthreads = 0L, 
                  const bool verbose = true) {
  
  const unsigned int N = actifs.size(), K = emplois.size();
  const std::size_t Ndata = xr_dist.size();

  // instantation de la fonction d'attraction.
  mode_attraction type_att;
  auto fct_attraction = type_att.create(parametres, attraction);
  
  Urban urb(jr_dist, p_dist, xr_dist, actifs, emplois, fuites);

   // Passage explicite en std::vector pour rendre les vecteurs thread safe (ts_)(nécessaire pour openmp dans la macro).
   //const std::vector<double> ts_emplois = as<std::vector<double>>(emplois);
   //const std::vector<double> ts_fuite = as<std::vector<double>>(fuites);
   //const std::vector<double> ts_actifs = as<std::vector<double>>(actifs);
   
  // const std::vector<double> ts_parametres = as< std::vector<double> >(parametres);
   
  // const std::vector<int> ts_jr_dist = as< std::vector<int> >(jr_dist);
  // const std::vector<int> ts_p_dist = as< std::vector<int> >(p_dist);
  // const std::vector<double> ts_xr_dist = as< std::vector<double> >(xr_dist);
   
  // const std::vector<double> ts_xr_odds = as< std::vector<double> >(xr_odds);
   
   std::vector<double> emplois_libres(emplois);
   std::vector<double> actifs_libres(actifs);
   
   // Initialisation de la matrice origines-destination.
   std::vector< std::vector<float> > liaisons(N, std::vector<float> (K));
   
   double tot_actifs_libres = std::accumulate(actifs_libres.begin(), actifs_libres.end(), 0.0);
   double tot_actifs = tot_actifs_libres, old_tot;
   
  
int nloop = 0;

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

#pragma omp parallel num_threads(ntr) 
{
  #pragma omp declare reduction(vsum : std::vector<double> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), \
  std::plus<double>())) initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

do { // DEBUT DES BOUCLES DO-WHILE
  
  #pragma omp single
  {
    nloop++;
    if (verbose == TRUE) REprintf("\nBoucle %i: ", nloop);
  }
  #pragma omp for schedule(static, 2)
  for (auto from = 0; from < N; ++from) {
    
    if (actifs_libres[from] < LIMITE_PRECISION_3) continue;
    
    Urban::Residents res(urb, from);
    unsigned int n_sites = res.col_dispo.size();
    std::vector<double> attirances(n_sites), cumul_att(n_sites);

    attirances = attractivite(res, fct_attraction);
    cumul_att = attract_cumul(res, attirances);

    // Attention : attirances stocke temporairement le résultat de la répartition des emplois.
    attirances = calc_repartition(res, attirances, cumul_att);

    for (auto k: res.col_dispo) {
      index = std::distance(res.col_dispo.begin(), k);
      liaisons[from][*k] += (float) attirances[index];
    }
  } // fin des boucles sur les from
  
  // Après la distribution simultanée, il faut récupérer les excédents et reconstruire des vecteurs
  // d'actifs encore libres et d'emplois encore libres.
  #pragma omp for 
  for (auto i = 0; i < N; ++i) {
    actifs_libres[i] = 0;
  }
  
  // Traitement des dépassements en colonnes.
  #pragma omp for reduction(vsum: actifs_libres)
  for (auto j = 0; j < K; ++j) {
    emplois_libres[j] = urb.emplois[j];
    for (auto i = 0; i < N; ++i) {
      emplois_libres[j] -= (double) liaisons[i][j];
    }
    if (emplois_libres[j] < 0) {
      double tx_depassement = urb.emplois[j] / (urb.emplois[j] - emplois_libres[j]);
      
      // Renvoie à domicile des actifs excédentaires à proportion de leur contribution à l'excédent.
      // CHOIX METHODO : le renvoi à domicile est à proportion du total des actifs en place, et non pas 
      // uniquement à proportion de la dernière vague d'arrivants.
      for (auto i = 0; i < N; ++i) {
        actifs_libres[i] += (double) liaisons[i][j] * (1 - tx_depassement);
        liaisons[i][j] = liaisons[i][j] * (float) tx_depassement; 
      }
      emplois_libres[j] = 0;
    }
  }
  
  #pragma omp single
  {   
    if (verbose == TRUE) 
      REprintf("%.0f actifs non occupés (soit %.1f %%)", tot_actifs_libres, 100 * tot_actifs_libres/tot_actifs);
    
    old_tot = tot_actifs_libres;
    tot_actifs_libres = 0;
  }
  #pragma omp for reduction(+: tot_actifs_libres) 
  for (auto i = 0; i < N; ++i) {
    actifs_libres[i] /= (1 - urb.fuite[i]); // Ne pas oublier de renvoyer aussi à domicile les fuyards des actifs renvoyés.
    tot_actifs_libres += actifs_libres[i];
  }
  
} while (tot_actifs_libres/tot_actifs > LIMITE_PRECISION_1 &&
         std::abs(tot_actifs_libres - old_tot)/tot_actifs > LIMITE_PRECISION_2 && 
         nloop < LIMITE_LOOP); // FIN DES BOUCLES DO-WHILE

} // fin clause omp parallel

// Mise en forme du résultat selon présence ou non de groupes et d'une cible.
if (group_from == R_NilValue || group_to == R_NilValue) {
   
  // Format de sortie
  std::vector<float> res_xr(Ndata);
  std::vector<int> res_i(Ndata);

  // sortie au format xr_dist.
  #pragma omp parallel for
  for (auto i = 0; i < N; ++i) {
    for (auto k = urb.p_dist[i]; k < urb.p_dist[i + 1]; ++k) {
      res_xr[k] = liaisons[i][ urb.jr_dist[k] ];
      res_i[k] = i;
    }
  }
  return List::create(_("i") = wrap(res_i), _("j") = jr_dist, _("flux") = wrap(res_xr));
} else {
  // sortie du résultat agrégé.
  const std::vector<int> ts_group_from = as< std::vector<int> >(group_from);
  const std::vector<int> ts_group_to = as< std::vector<int> >(group_to);
  int Nref = ts_group_from.size(), Kref = ts_group_to.size();
  
  NumericVector res_i(Nref * Kref), res_j(Nref * Kref);
  std::vector<double> flux(Nref * Kref);

  #pragma omp declare reduction(vsum : std::vector<double> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), \
  std::plus<double>())) initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
  #pragma omp parallel for reduction(vsum: flux) collapse(2)
  for (auto i = 0; i < N; ++i) {
    for (auto j = 0; j < K; ++j) {
      flux[ ts_group_from[i] * Kref + ts_group_to[j] ] += liaisons[i][j];
    }
  }
  
  for (auto i = 0; i < Nref; ++i) {
    for (auto j = 0; j < Kref; ++j) {
      res_i[i * Kref + j] = i;
      res_j[i * Kref + j] = j;
    }
  }
  
  if (cible == R_NilValue) {
    return List::create(_("i") = res_i, _("j") = res_j, _("flux") = wrap(flux));
  } else {
    std::vector<double> p_cible = as< std::vector<double> >(cible);
    std::vector<double> p_flux(flux);
    // Calcul du KL
    double tot_flux = std::accumulate(flux.begin(), flux.end(), 0.0);
    double tot_cible = std::accumulate(p_cible.begin(), p_cible.end(), 0.0);
    
    if (tot_flux == 0) stop("Les flux groupés sont tous nuls.");
    
    double kl = 0;
    
    for (auto k = 0; k < Nref * Kref; ++k) {
      p_cible[k] /= tot_cible;
      p_flux[k] /= tot_flux;
      if (p_flux[k] < PLANCHER_KL) p_flux[k] = PLANCHER_KL;
      kl += p_flux[k] * (log(p_flux[k]) - log(p_cible[k]));
    }
    
    return List::create(_("i") = res_i, _("j") = res_j, _("flux") = wrap(flux), _("kl") = kl);
  }
}
}
 
 