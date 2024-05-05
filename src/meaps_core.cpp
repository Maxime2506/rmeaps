#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <iterator>

#include "fcts_penal.h"

std::vector< std::vector<double> > meaps_core(const std::vector<int> jr_dist, 
                                              const std::vector<int> p_dist, 
                                              const std::vector<double> xr_dist, 
                                              std::vector<double> emplois,
                                              const std::vector<double> actifs, 
                                              std::vector<double> fuite, 
                                              const std::vector<double> parametres,
                                              const std::vector<double> xr_odds,
                                              const std::string attraction,
                                              const int nthreads, 
                                              const bool verbose) {
  const int _LIMITE_LOOP = 200; // condition d'arrêt pour les boucles lors de la distribution des résidents vers des emplois.
  const double LIMITE_PRECISION_1 = 1e-3; // condition d'arrêt sur le pourcentage de résidents non classés restants.
  const double LIMITE_PRECISION_2 = 1e-4; // condition d'arrêt sur la vitesse de reclassement des résidents non occupés.
  
  const std::size_t N = actifs.size(), K = emplois.size();
   
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
   
   std::vector<double> emplois_libres(emplois);
   std::vector<double> actifs_libres(actifs);
   
   // Initialisation du résultat.
   std::vector< std::vector<double> > liaisons(N, std::vector<double> (K));

   double tot_actifs_libres = std::accumulate(actifs_libres.begin(), actifs_libres.end(), 0.0);
   double tot_actifs = tot_actifs_libres, old_tot;
   
   int nloop = 0;
   
#pragma omp parallel num_threads(ntr) 
{
#pragma omp declare reduction(vsum : std::vector<double> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), \
  std::plus<double>())) initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
  do {
    
#pragma omp single
{
  nloop++;
  if (verbose == TRUE) REprintf("\nBoucle %i: ", nloop);
}
#pragma omp for schedule(static, 2)
for (std::size_t from = 0; from < N; ++from) {
  
  // Inner index issu de la matrice sparse pour la ligne from en cours.
  std::size_t debut = p_dist[from], fin = p_dist[from + 1L];
  std::size_t k_valid = fin - debut;
  
  // Calcul de l'attractivité d'un site modulée par la fonction d'attraction retenue.
  std::vector<double> attirances(k_valid), repartition(k_valid);
  
  for (std::size_t k = 0; k < k_valid; ++k) {
    attirances[k] = emplois_libres[ jr_dist[ debut + k] ];
  }
  
  if (attraction == "marche") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= marche(xr_dist[debut + k], parametres[0], parametres[1]);
    }}
  
  if (attraction == "marche_liss") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= marche_liss(xr_dist[debut + k], parametres[0], parametres[1]);
    }}
  
  if (attraction == "double_marche_liss") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= marche_liss(xr_dist[debut + k], parametres[0], parametres[1], parametres[2], parametres[3]);
    }}
  
  if (attraction == "decay") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= decay(xr_dist[debut + k], parametres[0], parametres[1]);
    }}
  
  if (attraction == "logistique") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= logistique(xr_dist[debut + k], parametres[0], parametres[1], parametres[2]);
    }}
  
  if (attraction == "odds") {
    for (std::size_t k = 0; k < k_valid; ++k) {
        attirances[k] *= xr_odds[debut + k];
      } 
    } 
  
  // Calcul de l'accessibilité pondérée par l'attraction.
  std::vector<double> accessibility(k_valid);
  double tot = 0;
  for (std::size_t k = 0; k < k_valid;) {
    auto pos = k + 1L;
    while (xr_dist[debut + k] == xr_dist[debut + pos] && pos < k_valid) ++pos;
    for (std::size_t ego = k; ego < pos; ++ego) {
      tot += attirances[ego];
    }
    for (std::size_t ego = k; ego < pos; ++ego) {
      accessibility[ego] = tot;
    }
    k = pos;
  }
  
  if (accessibility[k_valid - 1L] <= 0) continue;
  
  // Calcul de l'absorption sur la ligne from considérée.
  double absorption = -log(fuite[from]) / accessibility[k_valid - 1L];
  
  // Calcul des actifs absorbés par sites.
  std::vector<double> jobtakers(k_valid + 1L, actifs_libres[from]);
  for(std::size_t k = 0L; k < k_valid; ++k) {
    jobtakers[k + 1L] *= exp(-absorption * accessibility[k]); // ceux qui dépassent le site k+1.
  }
  for(std::size_t k = 0L; k < k_valid; ++k) {
    jobtakers[k] -= jobtakers[k + 1L];
  }
  // Répartition des jobtakers en traitant les cas à distances égales.
  for (std::size_t k = 0; k < k_valid;) {
    double tot_attirances = 0, tot_jobtakers = 0;
    auto pos = k + 1L;
    while (xr_dist[debut + k] == xr_dist[debut + pos] && pos < k_valid) ++pos;
    for (std::size_t ego = k; ego < pos; ++ego) {
      tot_attirances += attirances[ego];
      tot_jobtakers += jobtakers[ego];
    }
    if (tot_attirances > 0) {
      for (std::size_t ego = k; ego < pos; ++ego) {
        liaisons[from][ jr_dist[debut + ego] ] += attirances[ego] / tot_attirances * tot_jobtakers;
      }
    }
    k = pos;
  }
} 

// Après la distribution simultanée, il faut récupérer les excédents et reconstruire des vecteurs
// d'actifs encore libres et d'emplois encore libres.
#pragma omp for 
for (std::size_t i = 0; i < N; ++i) {
  actifs_libres[i] = 0;
}

// Traitement des dépassements en colonnes.
#pragma omp for reduction(vsum: actifs_libres)
for (std::size_t j = 0; j < K; ++j) {
  emplois_libres[j] = emplois[j];
  for (std::size_t i = 0; i < N; ++i) {
    emplois_libres[j] -= liaisons[i][j];
  }
  if (emplois_libres[j] < 0) {
    double tx_depassement = emplois[j] / (emplois[j] - emplois_libres[j]);
    
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
#pragma omp single
{   
  if (verbose == TRUE) 
    REprintf("%f actifs non occupés (soit %f %%)", round(tot_actifs_libres), round(1000 * tot_actifs_libres/tot_actifs)/10);
  
  old_tot = tot_actifs_libres;
  tot_actifs_libres = 0;
}
#pragma omp for reduction(+: tot_actifs_libres) 
for (std::size_t i = 0; i < N; ++i) {
  actifs_libres[i] /= (1 - fuite[i]); // Ne pas oublier de renvoyer aussi à domicile les fuyards des actifs renvoyés.
  tot_actifs_libres += actifs_libres[i];
}

  } while (tot_actifs_libres/tot_actifs > LIMITE_PRECISION_1 && 
    std::abs(tot_actifs_libres - old_tot)/tot_actifs > LIMITE_PRECISION_2 && 
    nloop < _LIMITE_LOOP);
} // fin clause omp parallel

return liaisons;
 }
 
 