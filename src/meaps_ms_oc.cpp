# ifdef _OPENMP
#include <omp.h>
# endif
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <algorithm>
#include <iterator>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
// #include "repartir_continu.h"
#include "classes_attraction.h"

// #include "fcts_penal.h"
using namespace Rcpp;

inline std::vector<double> _one_distrib_continu(
    const double entrants, 
    const double fuite,
    const std::vector<double>& attractivite,
    const std::vector<double>& distances,
    const std::vector<double>& placeslibres) {
  
  // Le cas sum(placeslibres) = 0 doit être géré en amont.
  std::size_t k_valid = placeslibres.size();
  std::vector<double> attraction(k_valid), accessibility(k_valid), placesprises(k_valid);
  
  // Calcul de l'attraction de chaque site.
  double total_places = 0;
  for (std::size_t k = 0; k < k_valid; ++k) {
    attraction[k] = placeslibres[k] * attractivite[k];
    total_places += placeslibres[k];
  }
  
  // Cas où il y a plus d'entrants que de places libres.
  if (total_places <= entrants) {
    return placeslibres;
  }
  
  // Calcul de l'accessibilité (corrigée de l'attraction).
  double tot = 0;
  for (std::size_t k = 0; k < k_valid;) {
    auto pos = k + 1L;
    while (distances[k] == distances[pos] && pos < k_valid) ++pos;
    for (std::size_t ego = k; ego < pos; ++ego) {
      tot += attraction[ego];
    }
    for (std::size_t ego = k; ego < pos; ++ego) {
      accessibility[ego] = tot;
    }
    k = pos;
  }
  
  // Cas (improbable) où l'accessibilité resterait nulle. En ce cas, pas de remplissage.
  if (accessibility[k_valid - 1L] == 0) {
    std::vector<double> zeros(k_valid);
    return zeros;
  }
  
  // Calcul de l'absorption.
  double absorption = -log(fuite) / accessibility[k_valid - 1L];
  
  // Calcul des actifs absorbés par sites (avant traitement des sites à distances égales).
  std::vector<double> jobtakers(k_valid + 1L);
  jobtakers[0L] = entrants;
  for(std::size_t k = 0L; k < k_valid; ++k) {
    jobtakers[k + 1L] = entrants * exp(-absorption * accessibility[k]); // ceux qui dépassent le site k+1.
  }
  for(std::size_t k = 0L; k < k_valid; ++k) {
    jobtakers[k] -= jobtakers[k + 1L];
  }
  
  // Répartition des jobtakers dans tous les cas.
  for (std::size_t k = 0; k < k_valid;) {
    double tot_attraction = 0;
    double tot_jobtakers = 0;
    auto pos = k + 1L;
    while (distances[k] == distances[pos] && pos < k_valid) ++pos;
    for (std::size_t ego = k; ego < pos; ++ego) {
      tot_attraction += attraction[ego];
      tot_jobtakers += jobtakers[ego];
    }
    if (tot_attraction > 0) {
      for (std::size_t ego = k; ego < pos; ++ego) {
        placesprises[ego] = attraction[ego] / tot_attraction * tot_jobtakers;
      }
    }
    k = pos;
  }
  
  return placesprises;
}

// Fonction de répartition des actifs entre les sites d'emplois selon l'attractivité du site et la fuite.
// Cette fonction gère le cas d'un dépassement de l'offre.
inline std::vector<double> _repartir_continu(
    const double actifs, 
    const double fuite,
    const std::vector<double>& attractivite,
    const std::vector<double>& distances,
    std::vector<double>& placeslibres) {
  
  std::size_t k_valid = placeslibres.size();
  double actifs_non_etablis = actifs, tot_placeslibres = 0.0;
  
  // Cas où il n'y a pas assez de places libres.
  for (std::size_t k = 0; k < k_valid; ++k) tot_placeslibres += placeslibres[k];
  if (actifs * (1 - fuite) >= tot_placeslibres) return placeslibres;
  
  std::vector<double> placesprises(k_valid, 0.0), placesrestantes(placeslibres), nouvellesprises(k_valid);
  
  do {
    nouvellesprises = _one_distrib_continu(actifs_non_etablis, fuite, attractivite, distances, placesrestantes);
    
    actifs_non_etablis = 0.0;
    for (std::size_t k = 0; k < k_valid; ++k) {
      placesprises[k] += nouvellesprises[k];
      if (placesprises[k] > placeslibres[k]) {
        actifs_non_etablis += placesprises[k] - placeslibres[k];
        placesprises[k] = placeslibres[k];
      }
      placesrestantes[k] = placeslibres[k] - placesprises[k];
    }
    // Attention : les actifs fuyant ne doivent pas être oubliés à côté de ceux ayant débordés. 
    // Il faut ajouter une part résiduelle de fuyards.
    actifs_non_etablis *= (1 + fuite);
  } while (actifs_non_etablis > 1e-9);// Des boucles qui n'apportent rien peuvent être effectuées si les conditions sont trop stricts. 
  
  return placesprises;
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
 
 // [[Rcpp::export]]
 
 List multishuf_oc_group_cpp(
     const IntegerVector jr_dist, 
     const IntegerVector p_dist, 
     const NumericVector xr_dist, 
     const NumericVector emplois,
     const NumericVector actifs, 
     const NumericVector fuites, 
     const IntegerMatrix shuf, 
     const IntegerVector group_from,
     const IntegerVector group_to,
     const NumericVector parametres,
     const NumericVector xr_odds,
     const Nullable<NumericVector> cible = R_NilValue,
     const std::string attraction = "constant",
     int nthreads = 0, bool verbose = true) {
   
   const std::size_t N = actifs.size(), Nboot = shuf.nrow(), Ns = shuf.ncol();
   const double PLANCHER_KL = 1e-6; // gestion de cases nulles dans le calcul de l'entropie relative (KL).
   // Instantation de la fonction d'attraction.
  const std::vector<double> param = as<std::vector<double> >(parametres);
  mode_attraction type_att;
  auto fct_attraction = type_att.create(param, attraction);

   auto Nref = *std::max_element(group_from.begin(), group_from.end());
   Nref = Nref + 1L;
   auto Kref = *std::max_element(group_to.begin(), group_to.end());
   Kref = Kref + 1L;
   int ntr = 1;
#ifdef _OPENMP
   ntr = nthreads;
   if (ntr == 0) {
     ntr = omp_get_max_threads();
   }
   if (ntr > omp_get_max_threads()) {
     ntr = omp_get_max_threads();
   }
   if(ntr>Nboot) ntr = Nboot;
   if (verbose == TRUE) REprintf("Nombre de threads = %i\n", ntr);
#endif
   
   // REMARQUE : les conversions de Numeric et IntegerVector en équivalent std::vector sont ici nécessaires car
   // les vecteurs initiaux de sont pas thread safe dans openmp.
   // Conversion en vecteur de vecteurs C++.
   std::vector< std::vector<int> > ishuf(Nboot, std::vector<int>(Ns));
   for (std::size_t i = 0; i < Nboot; ++i) {
     for (std::size_t j = 0; j < Ns; ++j) {
       ishuf[i][j] = shuf(i, j);
       if (ishuf[i][j] >= N) {
         stop("Les rangs de shufs vont au-delà de la taille du vecteur actifs.");
       }
     }
   }
   
   std::vector<double> xr_lor = as< std::vector<double> >(xr_odds);
   
   std::vector<double> emploisinitial = as<std::vector<double>>(emplois);
   std::vector<double> fcpp = as<std::vector<double>>(fuites);
   std::vector<double> actifscpp = as<std::vector<double>>(actifs);
   
   std::vector<int> _row_group = as< std::vector<int> >(group_from);
   std::vector<int> _col_group = as< std::vector<int> >(group_to);
   
   const std::vector<int> _jr_dist = as< std::vector<int> >(jr_dist);
   const std::vector<int> _p_dist = as< std::vector<int> >(p_dist);
   const std::vector<double> _xr_dist = as< std::vector<double> >(xr_dist);
   
   std::vector<double> resultat(Nref*Kref, 0);
   
   // Le vecteur shuf peut être plus long que le nombre de lignes de rkdist
   // s'il fait repasser plusieurs fois la même ligne d'actifs. Dans ce cas, on
   // compte la fréquence de passage de chaque ligne et l'on divise le poids de
   // la ligne par cette fréquence.
   std::vector<int> freq_actifs(N, 0L);
   for (auto i : ishuf[0]) {
     freq_actifs[i]++;
   }
   
   // Initialisation du résultat.
   // Lancement du bootstrap.
   Progress p(Nboot * Ns, verbose);
   std::size_t NNboot = floor( Nboot / ntr);
   if(NNboot == 0) NNboot = 1;
   // Un vecteur représentant la matrice des flux groupés.
   std::vector< std::vector<double> > liaisons(ntr, std::vector<double> (Nref * Kref));
   for(size_t Iboot = 0; Iboot < ntr; ++Iboot) {
     std::fill( liaisons[Iboot].begin(), liaisons[Iboot].end(), 0 );
   }
   
#pragma omp parallel num_threads(ntr)
{
#pragma omp for schedule(static)
  for(size_t Iboot = 0; Iboot < ntr; ++Iboot) {
    size_t max_i = (Iboot+1)*NNboot;
    if(Iboot==ntr-1) max_i = Nboot;
    size_t min_i = Iboot*NNboot;
    for (size_t iboot = min_i; iboot < max_i; ++iboot) {
      // Initialisation de l'ordre de départ des actifs et de l'emploi
      // disponible au début.
      
      std::vector<double> emp(emploisinitial);  // deep copy.
      
      for (auto from : ishuf[iboot]) {
        
        // Construction de l'accessibilité dite pénalisée.
        std::size_t debut = _p_dist[from], fin = _p_dist[from + 1L];
        std::size_t k_valid = fin - debut;
        std::vector<double> facteur_attraction(k_valid, 1.0), 
        emplois_libres(k_valid),
        repartition(k_valid, 0);
        
        for (std::size_t k = 0; k < k_valid; ++k) {
          facteur_attraction[k]  = (*fct_attraction)(_xr_dist[debut + k]);
        }
        
        for (std::size_t k = 0; k < k_valid; ++k) {
          emplois_libres[k] = emp[ _jr_dist[ debut + k] ];
        }
        
        double actifspartant = actifscpp[from] / freq_actifs[from];
        
        std::vector<double> dist(_xr_dist.begin() + debut, _xr_dist.begin() + fin);
        repartition = _repartir_continu(actifspartant, fcpp[from], facteur_attraction, dist, emplois_libres);
        
        // Impact sur l'emploi disponible total et sommation sur les emplois pris.
        std::size_t curseur_ligne = _row_group[from] * Kref;
        
        for (std::size_t k = 0; k < k_valid; ++k) {
          emp[ _jr_dist[debut + k] ] -= repartition[k];
          liaisons[Iboot][ curseur_ligne + _col_group[_jr_dist[debut + k]] ] += repartition[k];
        }
        // check interrupt & progress
        p.increment();
        if (from % 100 == 0) Progress::check_abort();
      } // shuf
    } // iboot
  } // Iboot
#pragma omp for
  for (std::size_t i = 0; i < Nref; ++i) {
    for (std::size_t j = 0; j < Kref; ++j) {
      for(size_t Iboot = 0; Iboot < ntr; ++Iboot) {
        resultat[i * Kref + j] += liaisons[Iboot][i * Kref + j ] / Nboot;
      }
    }
  }
} // omp

NumericVector out(Nref*Kref);
IntegerVector res_i(Nref*Kref), res_j(Nref*Kref);
for (std::size_t i = 0; i < Nref; ++i) {
  for (std::size_t j = 0; j < Kref; ++j) {
    for(size_t Iboot = 0; Iboot < ntr; ++Iboot) {
      out(i*Kref + j) += liaisons[Iboot][i * Kref + j ] / Nboot;
    }
    res_i(i*Kref + j) = i;
    res_j(i*Kref + j) = j;
  }
}

if (cible.isNull()) {
  return List::create(_("i") = res_i, _("j") = res_j, _("flux") = out);
} 

std::vector<double> p_cible = as< std::vector<double> >(cible);
std::vector<double> p_flux= as< std::vector<double> >(out);
if(p_cible.size() != p_flux.size())
  stop("La cible n'a pas la bonne longueur");
// Calcul du KL
double tot_flux = std::accumulate(out.begin(), out.end(), 0.0);
double tot_cible = std::accumulate(p_cible.begin(), p_cible.end(), 0.0);

if (tot_flux == 0) stop("Les flux groupés sont tous nuls");

double kl = 0;
double kl_term;
for (auto k = 0; k < Nref * Kref; ++k) {
  if (p_flux[k] == 0) {
    continue; 
  }
  p_cible[k] /= tot_cible;
  p_flux[k] /= tot_flux;
  if (p_cible[k] < PLANCHER_KL) p_cible[k] = PLANCHER_KL;
  kl_term = p_flux[k] * (log(p_flux[k]) - log(p_cible[k]));
  kl += kl_term;
}

return List::create(_("i") = res_i, _("j") = res_j, _("flux") = out, _("kl") = kl);
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
  for (auto i = 0; i < cible.nrow(); ++i){
    for (auto j = 0; j < cible.ncol(); ++j) {
      if (cible(i,j) != 0 && estim(i,j) != 0) tot += estim(i,j) * (log(estim(i,j)) - log(cible(i,j)));
    }
  }
  
  return tot;
}

//'
 // [[Rcpp::export]]
 IntegerVector max_threads() {
   return(wrap(omp_get_max_threads()));
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
 
 // [[Rcpp::export]]
 List multishuf_oc_cpp(
     const IntegerVector jr_dist, 
     const IntegerVector p_dist, 
     const NumericVector xr_dist, 
     const NumericVector emplois,
     const NumericVector actifs, 
     const NumericVector fuites, 
     const IntegerMatrix shuf, 
     const NumericVector parametres,
     const NumericVector xr_odds,
     const std::string attraction = "constant",
     int nthreads = 0, bool verbose = true) {
   
   const std::size_t N = actifs.size(), Nx  = xr_dist.size(),
     Nboot = shuf.nrow(), Ns = shuf.ncol();
    // Instantation de la fonction d'attraction.
  const std::vector<double> param = as<std::vector<double> >(parametres);
  mode_attraction type_att;
  auto fct_attraction = type_att.create(param, attraction);

   int ntr = 1;
#ifdef _OPENMP
   ntr = nthreads;
   if (ntr == 0) {
     ntr = omp_get_max_threads();
   }
   if (ntr > omp_get_max_threads()) {
     ntr = omp_get_max_threads();
   }
   if(ntr>Nboot) ntr = Nboot;
   if (verbose == TRUE) REprintf("Nombre de threads = %i\n", ntr);
#endif
   
   // REMARQUE : les conversions de Numeric et IntegerVector en équivalent std::vector sont ici nécessaires car
   // les vecteurs initiaux de sont pas thread safe dans openmp.
   // Conversion en vecteur de vecteurs C++.
   std::vector< std::vector<int> > ishuf(Nboot, std::vector<int>(Ns));
   for (std::size_t i = 0; i < Nboot; ++i) {
     for (std::size_t j = 0; j < Ns; ++j) {
       ishuf[i][j] = shuf(i, j);
       if (ishuf[i][j] >= N) {
         stop("Les rangs de shufs vont au-delà de la taille du vecteur actifs.");
       }
     }
   }
   
   std::vector<double> xr_lor = as< std::vector<double> >(xr_odds);
   
   std::vector<double> emploisinitial = as<std::vector<double>>(emplois);
   std::vector<double> fcpp = as<std::vector<double>>(fuites);
   std::vector<double> actifscpp = as<std::vector<double>>(actifs);
      
   const std::vector<int> _jr_dist = as< std::vector<int> >(jr_dist);
   const std::vector<int> _p_dist = as< std::vector<int> >(p_dist);
   const std::vector<double> _xr_dist = as< std::vector<double> >(xr_dist);
   const std::vector<double> _xr_odds = as< std::vector<double> >(xr_odds);
   
   // 
   // Le vecteur shuf peut être plus long que le nombre de lignes de rkdist
   // s'il fait repasser plusieurs fois la même ligne d'actifs. Dans ce cas, on
   // compte la fréquence de passage de chaque ligne et l'on divise le poids de
   // la ligne par cette fréquence.
   std::vector<int> freq_actifs(N, 0L);
   for (auto i : ishuf[0]) {
     freq_actifs[i]++;
   }
   
   // Initialisation du résultat.
   // Lancement du bootstrap.
   Progress p(Nboot*Ns, verbose);
   
   std::size_t NNboot = floor( Nboot / ntr );
   if(NNboot == 0) NNboot = 1;
   // Un vecteur représentant la matrice des flux groupés.
   std::vector< std::vector<double> > liaisons(ntr, std::vector<double> (Nx));
   for(size_t Iboot = 0; Iboot < ntr; ++Iboot) {
     std::fill( liaisons[Iboot].begin(), liaisons[Iboot].end(), 0 );
   }
   
#pragma omp parallel num_threads(ntr) 
{
#pragma omp for schedule(static)
  for(size_t Iboot = 0; Iboot < ntr; ++Iboot) {
    size_t max_i = (Iboot+1)*NNboot;
    if(Iboot==ntr-1) max_i = Nboot;
    size_t min_i = Iboot*NNboot;
    if(min_i<max_i) {
      for (size_t iboot = min_i; iboot < max_i; ++iboot) {
        // Initialisation de l'ordre de départ des actifs et de l'emploi
        // disponible au début.
        
        
        std::vector<double> emp(emploisinitial);  // deep copy.
        
        for (auto from : ishuf[iboot]) {
          // check interrupt & progress
          p.increment();
          if(from % 100==0)
            Progress::check_abort();
          
          // Construction de l'accessibilité dite pénalisée.
          std::size_t debut = _p_dist[from], fin = _p_dist[from + 1L];
          std::size_t k_valid = fin - debut;
          std::vector<double> facteur_attraction(k_valid, 1.0), 
          emplois_libres(k_valid),
          repartition(k_valid, 0);
          
          for (std::size_t k = 0; k < k_valid; ++k) {
          facteur_attraction[k]  = (*fct_attraction)(_xr_dist[debut + k]);
          }
          
          for (std::size_t k = 0; k < k_valid; ++k) {
            emplois_libres[k] = emp[ _jr_dist[ debut + k] ];
          }
          
          double actifspartant = actifscpp[from] / freq_actifs[from];
          
          std::vector<double> dist(_xr_dist.begin() + debut, _xr_dist.begin() + fin);
          repartition = _repartir_continu(actifspartant, fcpp[from], facteur_attraction, dist, emplois_libres);
          
          for (std::size_t k = 0; k < k_valid; ++k) {
            emp[ _jr_dist[debut + k] ] -= repartition[k];
            liaisons[Iboot][ debut + k ] += repartition[k];
          }
        } // if min_i<max_i
      } // shuf
    } // iboot
  } // Iboot
} // omp parallel

NumericVector resultat(Nx);

for (std::size_t i = 0; i < Nx; ++i) {
  for(size_t Iboot = 0; Iboot < ntr; ++Iboot) {
    resultat[i] += liaisons[Iboot][i] / Nboot;
  }
}

// IntegerVector res_i(Nx);
// std::size_t ii = 0;
// for (std::size_t i = 0; i < N; ++i) {
//   for (std::size_t k = p_dist[i]; k < p_dist[i + 1L]; ++k) {
//     res_i[ii] = i;
//     ii++;
//   }
// }
return List::create(_("flux") = resultat);
// return List::create(_("i") = res_i,
//                     _("j") = jr_dist,
//                     _("flux") = resultat);
 }
