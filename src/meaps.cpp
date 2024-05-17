#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <iterator>

#include "fcts_penal.h"

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
 //' @param nthreads Nombre de threads pour OpenMP. Default : 0 = choix auto. 
 //' @param verbose Défaut = true. 
 //' On doit gérer normalisation et fuite_min auparavant.
 //'
 //' @return renvoie les flux au format triplet.
 // [[Rcpp::export]]
 DataFrame meaps_all_in(const IntegerVector jr_dist, 
                        const IntegerVector p_dist, 
                        const NumericVector xr_dist, 
                        NumericVector emplois,
                        const NumericVector actifs, 
                        NumericVector fuites, 
                        const NumericVector parametres,
                        const NumericVector xr_odds,
                        const std::string attraction = "constant",
                        const int nthreads = 0, 
                        const bool verbose = true) {
   
   const std::size_t N = actifs.size(), K = emplois.size(), Ndata = xr_dist.size();
   
   const int LIMITE_LOOP = 200; // condition d'arrêt pour les boucles lors de la distribution des résidents vers des emplois.
   const double LIMITE_PRECISION_1 = 1e-4; // condition d'arrêt sur le pourcentage de résidents non classés restants.
   const double LIMITE_PRECISION_2 = 1e-4; // condition d'arrêt sur la vitesse de reclassement des résidents non occupés.
   
   // Passage explicite en std::vector pour rendre les vecteurs thread safe (ts_)(nécessaire pour openmp dans la macro).
   std::vector<double> ts_emplois = as<std::vector<double>>(emplois);
   std::vector<double> ts_fuite = as<std::vector<double>>(fuites);
   const std::vector<double> ts_actifs = as<std::vector<double>>(actifs);
   
   const std::vector<double> ts_parametres = as< std::vector<double> >(parametres);
   
   const std::vector<int> ts_jr_dist = as< std::vector<int> >(jr_dist);
   const std::vector<int> ts_p_dist = as< std::vector<int> >(p_dist);
   const std::vector<double> ts_xr_dist = as< std::vector<double> >(xr_dist);
   
   const std::vector<double> ts_xr_odds = as< std::vector<double> >(xr_odds);
   
   std::vector<double> emplois_libres(ts_emplois);
   std::vector<double> actifs_libres(ts_actifs);
   
   // Initialisation de la matrice origines-destination.
   std::vector< std::vector<double> > liaisons(N, std::vector<double> (K));
   
   // Format de sortie
   std::vector<double> res_xr(Ndata);
   std::vector<int> res_i(Ndata);
   
   double tot_actifs_libres = std::accumulate(actifs_libres.begin(), actifs_libres.end(), 0.0);
   double tot_actifs = tot_actifs_libres, old_tot;
   
   const double LIMITE_PRECISION_3 = LIMITE_PRECISION_1 / tot_actifs; // seuil en deçà duquel une ligne d'actifs est supposée vides (tous actifs occupés). 
   
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
  // Inner index issu de la matrice sparse pour la ligne from en cours.
  std::size_t debut = ts_p_dist[from], fin = ts_p_dist[from + 1L];
  std::size_t k_valid = fin - debut;
  
  // Calcul de l'attractivité d'un site modulée par la fonction d'attraction retenue.
  std::vector<double> attirances(k_valid), repartition(k_valid);
  
  for (std::size_t k = 0; k < k_valid; ++k) {
    attirances[k] = emplois_libres[ ts_jr_dist[ debut + k] ];
  }
  
  if (attraction == "marche") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= marche(ts_xr_dist[debut + k], ts_parametres[0], ts_parametres[1]);
    }}
  
  if (attraction == "marche_liss") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= marche_liss(xr_dist[debut + k], ts_parametres[0], ts_parametres[1]);
    }}
  
  if (attraction == "double_marche_liss") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= double_marche_liss(xr_dist[debut + k], ts_parametres[0], ts_parametres[1], ts_parametres[2], ts_parametres[3]);
    }}
  
  if (attraction == "decay") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= decay(ts_xr_dist[debut + k], ts_parametres[0], ts_parametres[1]);
    }}
  
  if (attraction == "logistique") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= logistique(ts_xr_dist[debut + k], ts_parametres[0], ts_parametres[1], ts_parametres[2]);
    }}
  
  if (attraction == "odds") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= ts_xr_odds[debut + k];
    } 
  } 
  
  // Calcul de l'accessibilité pondérée par l'attraction.
  std::vector<double> accessibility(k_valid);
  double tot = 0;
  for (std::size_t k = 0; k < k_valid;) {
    auto pos = k + 1L;
    while (ts_xr_dist[debut + k] == ts_xr_dist[debut + pos] && pos < k_valid) ++pos;
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
  double absorption = -log(ts_fuite[from]) / accessibility[k_valid - 1L];
  
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
    while (ts_xr_dist[debut + k] == ts_xr_dist[debut + pos] && pos < k_valid) ++pos;
    for (std::size_t ego = k; ego < pos; ++ego) {
      tot_attirances += attirances[ego];
      tot_jobtakers += jobtakers[ego];
    }
    if (tot_attirances > 0) {
      for (std::size_t ego = k; ego < pos; ++ego) {
        liaisons[from][ ts_jr_dist[debut + ego] ] += attirances[ego] / tot_attirances * tot_jobtakers;
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
  emplois_libres[j] = ts_emplois[j];
  for (std::size_t i = 0; i < N; ++i) {
    emplois_libres[j] -= liaisons[i][j];
  }
  if (emplois_libres[j] < 0) {
    double tx_depassement = ts_emplois[j] / (ts_emplois[j] - emplois_libres[j]);
    
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
    REprintf("%.0f actifs non occupés (soit %.1f %%)", tot_actifs_libres, 100 * tot_actifs_libres/tot_actifs);
  
  old_tot = tot_actifs_libres;
  tot_actifs_libres = 0;
}
#pragma omp for reduction(+: tot_actifs_libres) 
for (std::size_t i = 0; i < N; ++i) {
  actifs_libres[i] /= (1 - ts_fuite[i]); // Ne pas oublier de renvoyer aussi à domicile les fuyards des actifs renvoyés.
  tot_actifs_libres += actifs_libres[i];
}

  } while (tot_actifs_libres/tot_actifs > LIMITE_PRECISION_1 &&
    std::abs(tot_actifs_libres - old_tot)/tot_actifs > LIMITE_PRECISION_2 && 
    nloop < LIMITE_LOOP); // FIN DES BOUCLES DO-WHILE
  

   // sortie au format xr_dist.
   
#pragma omp parallel for
   for (std::size_t i = 0; i < N; ++i) {
     for (std::size_t k = ts_p_dist[i]; k < ts_p_dist[i + 1L]; ++k) {
       res_xr[k] = liaisons[i][ ts_jr_dist[k] ];
       res_i[k] = i;
     }
   }
} // fin clause omp parallel
  if(verbose == TRUE) REprintf("\n");

   return DataFrame::create(_("i") = wrap(res_i), _("j") = jr_dist, _("flux") = wrap(res_xr));
 }
 
 
 //' La fonction meaps pour OPTIM qui distribue tous les actifs en même temps. En entrée, la matrice des distances (et si besoin des odds)
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
 NumericVector all_in_optim(const IntegerVector jr_dist, 
                            const IntegerVector p_dist, 
                            const NumericVector xr_dist, 
                            const IntegerVector group_from,
                            const IntegerVector group_to,
                            NumericVector emplois,
                            const NumericVector actifs, 
                            NumericVector fuites, 
                            const NumericVector parametres,
                            const NumericVector xr_odds,
                            const std::string attraction = "constant",
                            const int nthreads = 0, 
                            const bool verbose = true) {
   
   const int LIMITE_LOOP = 200; // condition d'arrêt pour les boucles lors de la distribution des résidents vers des emplois.
   const double LIMITE_PRECISION_1 = 1e-3; // condition d'arrêt sur le pourcentage de résidents non classés restants.
   const double LIMITE_PRECISION_2 = 1e-4; // condition d'arrêt sur la vitesse de reclassement des résidents non occupés.
   
   const std::size_t N = actifs.size(), K = emplois.size();
   auto Nref = 1L + *std::max_element(group_from.begin(), group_from.end());
   auto Kref = 1L + *std::max_element(group_to.begin(), group_to.end());

   // Passage explicite en std::vector pour rendre les vecteurs thread safe (ts_)(nécessaire pour openmp dans la macro).
   std::vector<double> ts_emplois = as<std::vector<double>>(emplois);
   std::vector<double> ts_fuite = as<std::vector<double>>(fuites);
   const std::vector<double> ts_actifs = as<std::vector<double>>(actifs);
   
   const std::vector<double> ts_parametres = as< std::vector<double> >(parametres);
   
   const std::vector<int> ts_jr_dist = as< std::vector<int> >(jr_dist);
   const std::vector<int> ts_p_dist = as< std::vector<int> >(p_dist);
   const std::vector<double> ts_xr_dist = as< std::vector<double> >(xr_dist);
   
   const std::vector<double> ts_xr_odds = as< std::vector<double> >(xr_odds);
   
   std::vector<double> emplois_libres(ts_emplois);
   std::vector<double> actifs_libres(ts_actifs);
   
   // Initialisation du résultat.
   std::vector< std::vector<double> > liaisons(N, std::vector<double> (K));
   
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
  do {
    
#pragma omp single
{
  nloop++;
  if (verbose == TRUE) REprintf("\nBoucle %i: ", nloop);
}
#pragma omp for schedule(static, 2)
for (std::size_t from = 0; from < N; ++from) {
  
  // Inner index issu de la matrice sparse pour la ligne from en cours.
  std::size_t debut = ts_p_dist[from], fin = ts_p_dist[from + 1L];
  std::size_t k_valid = fin - debut;
  
  // Calcul de l'attractivité d'un site modulée par la fonction d'attraction retenue.
  std::vector<double> attirances(k_valid), repartition(k_valid);
  
  for (std::size_t k = 0; k < k_valid; ++k) {
    attirances[k] = emplois_libres[ ts_jr_dist[ debut + k] ];
  }
  
  if (attraction == "marche") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= marche(ts_xr_dist[debut + k], ts_parametres[0], ts_parametres[1]);
    }
  }
  
  if (attraction == "marche_liss") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= marche_liss(xr_dist[debut + k], ts_parametres[0], ts_parametres[1]);
    }
  }
  
  if (attraction == "double_marche_liss") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= double_marche_liss(xr_dist[debut + k], ts_parametres[0], ts_parametres[1], ts_parametres[2], ts_parametres[3]);
    }
  }
  
  if (attraction == "decay") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= decay(ts_xr_dist[debut + k], ts_parametres[0], ts_parametres[1]);
    }
  }
  
  if (attraction == "logistique") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= logistique(ts_xr_dist[debut + k], ts_parametres[0], ts_parametres[1], ts_parametres[2]);
    }
  }
  
  if (attraction == "odds") {
    for (std::size_t k = 0; k < k_valid; ++k) {
      attirances[k] *= ts_xr_odds[debut + k];
    } 
  }
  
  // Calcul de l'accessibilité pondérée par l'attraction.
  std::vector<double> accessibility(k_valid);
  double tot = 0;
  for (std::size_t k = 0; k < k_valid;) {
    auto pos = k + 1L;
    while (ts_xr_dist[debut + k] == ts_xr_dist[debut + pos] && pos < k_valid) ++pos;
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
  double absorption = -log(ts_fuite[from]) / accessibility[k_valid - 1L];
  
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
    while (ts_xr_dist[debut + k] == ts_xr_dist[debut + pos] && pos < k_valid) ++pos;
    for (std::size_t ego = k; ego < pos; ++ego) {
      tot_attirances += attirances[ego];
      tot_jobtakers += jobtakers[ego];
    }
    if (tot_attirances > 0) {
      for (std::size_t ego = k; ego < pos; ++ego) {
        liaisons[from][ ts_jr_dist[debut + ego] ] += attirances[ego] / tot_attirances * tot_jobtakers;
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
  emplois_libres[j] = ts_emplois[j];
  for (std::size_t i = 0; i < N; ++i) {
    emplois_libres[j] -= liaisons[i][j];
  }
  if (emplois_libres[j] < 0) {
    double tx_depassement = ts_emplois[j] / (ts_emplois[j] - emplois_libres[j]);
    
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
    REprintf("%.0f actifs non occupés (soit %.1f %%)", tot_actifs_libres, 100 * tot_actifs_libres/tot_actifs);
  
  old_tot = tot_actifs_libres;
  tot_actifs_libres = 0;
}
#pragma omp for reduction(+: tot_actifs_libres) 
for (std::size_t i = 0; i < N; ++i) {
  actifs_libres[i] /= (1 - ts_fuite[i]); // Ne pas oublier de renvoyer aussi à domicile les fuyards des actifs renvoyés.
  tot_actifs_libres += actifs_libres[i];
}

  } while (tot_actifs_libres/tot_actifs > LIMITE_PRECISION_1 &&
    std::abs(tot_actifs_libres - old_tot)/tot_actifs > LIMITE_PRECISION_2 && 
    nloop < LIMITE_LOOP);
} // fin clause omp parallel

// sortie du résultat agrégé.
NumericVector resultat(Nref * Kref);
for (std::size_t i = 0; i < N; ++i) {
  for (std::size_t j = 0; j < K; ++j) {
    resultat(group_from[i] * Kref + group_to[j]) += liaisons[i][j];
  }
}

return resultat;
}
 


