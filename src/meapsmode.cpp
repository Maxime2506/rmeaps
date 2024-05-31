#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>

#include <algorithm>
#include <array>
#include <memory>

#include "classes.h"
#include "classes_attraction.h"
#include "constants.h"

using namespace Rcpp;

//' La fonction meaps qui distribue tous les actifs en même temps. En entrée, la matrice des distances (et si besoin des
//odds) ' doit être définie sous forme des inner et outer index d'une matrice sparse en ligne. Ceci revient à classer un
//data.frame ' avec les colonnes i, j et dist, d'abord par i, puis dist (=xr), puis j (=jr). ' @param jr_dist Le vecteur
//des indices des colonnes non vides. ' @param p_dist Le vecteur du nombres de valeurs non nulles sur chacune des
//lignes. ' @param xr_dist Le vecteur des valeurs dans l'ordre de jr_dist. ' @param emplois Le vecteur des emplois
//disponibles sur chacun des sites j (= marge des colonnes). ' @param actifs Le vecteur des actifs partant de chacune
//des lignes visées par shuf. Le vecteur doit faire la même longueur que shuf. ' @param fuite Le vecteur de la
//probabilité de fuite des actifs hors de la zone d'étude. ' @param parametres Un vecteur avec les paramètres
//nécessaires selon la fonction d'attraction retenue; ' @param attraction Choix de la fonction d'attraction des
//différents sites, appliquée à l'accessibilité. ' Par défaut, "constant" où aucun site n'a plus d'attrait qu'un autre.
//' "marche" où l'attrait vaut 1 jusqu'à une certaine distance (param 1) puis moins (param 2). f(x) = 1 si x < p1, = p2
//si x > p1. ' "logistique" où l'attrait décroît selon une fonction logistique avec une distance de bascule (param 1),
//une vitesse de bascule (param 2) ' et un seuil (param p). Si h(x) = exp( (x-p1)/p2), f(x) = p3 + h(x) / (1 + h(x)). '
//"odds" où chaque flux (from, to) se voit attribuer un odds. ' @param param est un vecteur avec dans l'ordre les
//valeurs des paramètres. ' @param group_from Le vecteur des regroupements des lignes d'actifs. ' @param group_to Le
//vecteur des regroupements des colonnes d'emplois. ' @param cible Le vecteur des flux groupés observés dans l'ordre des
//index c++ de (group_from, group_to). ' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto. '
//@param verbose Défaut = true. ' On doit gérer normalisation et fuite_min auparavant.
//'
//' @return renvoie les flux au format triplet.
// [[Rcpp::export]]
List meapsmode(const IntegerVector jr_dist, const IntegerVector p_dist, const NumericVector xr_dist,
                const NumericVector emplois, const NumericVector actifs, const NumericVector fuites,
                const IntegerVector j_mode, const NumericVector parametres, const std::string attraction = "constant",
                const Nullable<IntegerVector> group_from = R_NilValue,
                const Nullable<IntegerVector> group_to = R_NilValue, const Nullable<NumericVector> cible = R_NilValue,
                const int nthreads = 0L, const bool verbose = true) {
  const int N = actifs.size(), K = j_mode.size(), Nsites = emplois.size();

  // Instantation de la fonction d'attraction.
  const std::vector<double> param = as<std::vector<double> >(parametres);
  mode_attraction type_att;
  auto fct_attraction = type_att.create(param, attraction);

  Urban urb(jr_dist, p_dist, xr_dist, actifs, emplois, fuites);

  // Correspondance j (toidINS) vers j_mode.
  std::vector< std::vector<int> > col_mode(Nsites);
  for (auto j = 0; j < K; ++j) {
    col_mode[ j_mode[j] ].push_back(j);
  }

  std::vector<double> emplois_libres(urb.emplois);
  std::vector<double> emplois_distribues(K); // les emplois distribuées sur l'ensemble des colonnes avec mode.
  std::vector<double> actifs_libres(urb.actifs);

  for (auto j = 0; j < K; ++j) emplois_distribues[j] = emplois_libres[ j_mode[j] ];

  // Initialisation de la matrice origines-destination.
  std::vector<std::vector<float> > liaisons(N, std::vector<float>(K));

  double tot_actifs_libres = std::accumulate(urb.actifs.begin(), urb.actifs.end(), 0.0);
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

#pragma omp declare reduction(vsum : std::vector<double> : std::transform(                     \
        omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

#pragma omp parallel num_threads(ntr)
  {       // DEBUT DE LA PARALLELISATION
    do {  // DEBUT DES BOUCLES DO-WHILE

#pragma omp single
      {
        nloop++;
        if (verbose == TRUE) REprintf("Boucle %i: ", nloop);
      }
      
#pragma omp for schedule(static, 2)
      for (auto from = 0; from < N; ++from) {
        if (actifs_libres[from] < LIMITE_PRECISION_3) continue;

        Urban::Residents res(urb, from);

        // Repérage des colonnes disponibles rangées selon la distance pour les résidents res.
        std::vector<int> col_dispo = res.map_col_dispo(emplois_distribues);
        int n_sites = col_dispo.size();

        // Evaluation de l'attirance de chaque site.
        std::vector<double> attirances = res.attractivite(emplois_distribues, col_dispo, fct_attraction);

        // Puis répartition sans prise en compte des limites de places.
        attirances = res.repartition_nolimit(col_dispo, attirances, actifs_libres[from]);

        //
        for (auto k = 0; k < n_sites; ++k) {
          liaisons[from][col_dispo[k]] += static_cast<float>(attirances[k]);
        }
      }  // fin des boucles sur les from

// Après la distribution simultanée, il faut récupérer les excédents et reconstruire des vecteurs
// d'actifs encore libres et d'emplois encore libres.
#pragma omp for
      for (auto i = 0; i < N; ++i) actifs_libres[i] = 0;

// Traitement des dépassements en colonnes.
#pragma omp for reduction(vsum : actifs_libres)
      for (auto j = 0; j < Nsites; ++j) {
        double tot = 0;
        for (auto m: col_mode[j]) {
          emplois_distribues[m] = urb.emplois[j];
          for (auto i = 0; i < N; ++i) emplois_distribues[m] -= liaisons[i][m];
          tot += emplois_distribues[m];
        }

        if (tot > urb.emplois[j]) {
          double tx = urb.emplois[j] / tot;
          for (auto m: col_mode[j]) {
            for (auto i = 0; i < N; ++i) {
              actifs_libres[i] += (double)liaisons[i][m] * (1 - tx);
              liaisons[i][m] *= tx;
            }
            emplois_distribues[m] = 0;
          }
        }
      }

#pragma omp single
      {
        if (verbose == TRUE)
          REprintf("%.0f actifs non occupés (soit %.1f %%)\n", tot_actifs_libres, 100 * tot_actifs_libres / tot_actifs);

        old_tot = tot_actifs_libres;
        tot_actifs_libres = 0;
      }
#pragma omp for reduction(+ : tot_actifs_libres)
      for (auto i = 0; i < N; ++i) {
        actifs_libres[i] /= (1 - urb.fuites[i]);  // Ne pas oublier de renvoyer aussi à domicile les fuyards des actifs renvoyés.
        tot_actifs_libres += actifs_libres[i];
      }

    } while (tot_actifs_libres / tot_actifs > LIMITE_PRECISION_1 &&
             std::abs(tot_actifs_libres - old_tot) / tot_actifs > LIMITE_PRECISION_2 &&
             nloop < LIMITE_LOOP);  // FIN DES BOUCLES DO-WHILE

  }  // fin clause omp parallel

  if (verbose == TRUE) {
    REprintf("conditions finales :\n");
    REprintf("1 - Proportion des actifs non occupés = %.3e (seuil : %.3e)\n", tot_actifs_libres / tot_actifs,
             LIMITE_PRECISION_1);
    REprintf("2 - Vitesse de convergence = %.3e (seuil : %.3e)\n", std::abs(tot_actifs_libres - old_tot) / tot_actifs,
             LIMITE_PRECISION_2);
    REprintf("3 - Nombre de boucles = %i (seuil : %i )", nloop, LIMITE_LOOP);
  }
  // Mise en forme du résultat selon présence ou non de groupes et d'une cible.
  if (group_from.isNull() || group_to.isNull()) {
    return urb.format_sortie(liaisons);  // Retour de la matrice détaillée au format triplet.

  } else {
    std::vector<int> gfrom = as<std::vector<int> >(group_from);
    std::vector<int> gto = as<std::vector<int> >(group_to);

    Urban::SubRegion districts(urb, gfrom, gto);

    if (cible.isNull()) {
      return districts.format_sortie(liaisons);  // Retour de la matrice agrégrée au format triplet.
    } else {
      std::vector<float> ref = as<std::vector<float> >(cible);
      return districts.format_sortie(liaisons, ref);  // Retour de la matrice agrégrée au format triplet + KL.
    }
  }
}
