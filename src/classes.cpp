#ifdef _OPENMP
#include <omp.h>
#endif
#include "classes.h"

#include <algorithm>
#include <iterator>
#include <vector>

#include "classes_attraction.h"
#include "constants.h"
#include "entropie.h"

// Constructeurs Urban
Urban::Urban(){};

Urban::Urban(Rcpp::IntegerVector jr_, Rcpp::IntegerVector p_, Rcpp::NumericVector xr_, Rcpp::NumericVector actifs_,
             Rcpp::NumericVector emplois_, Rcpp::NumericVector fuites_)
    : jr(Rcpp::as<std::vector<int> >(jr_)),
      p(Rcpp::as<std::vector<int> >(p_)),
      xr(Rcpp::as<std::vector<double> >(xr_)),
      actifs(Rcpp::as<std::vector<double> >(actifs_)),
      emplois(Rcpp::as<std::vector<double> >(emplois_)),
      fuites(Rcpp::as<std::vector<double> >(fuites_)) {};

Urban::Urban(const std::vector<int> jr_, const std::vector<int> p_, const std::vector<double> xr_,
             const std::vector<double> actifs_, const std::vector<double> emplois_, const std::vector<double> fuites_)
    : jr(jr_), p(p_), xr(xr_), actifs(actifs_), emplois(emplois_), fuites(fuites_) {};

// Sortie : produit le triplet de la matrice détaillée des flux estimées par meaps.
Rcpp::List Urban::format_sortie(std::vector<std::vector<float> >& liens) {
  std::vector<float> res_xr(xr.size());
  std::vector<int> res_i(xr.size());

#pragma omp parallel for
  for (auto i = 0; i < n_from(); ++i) {
    for (auto k = p[i]; k < p[i + 1]; ++k) {
      res_xr[k] = liens[i][jr[k]];
      res_i[k] = i;
    }
  }
  return Rcpp::List::create(Rcpp::Named("i") = Rcpp::wrap(res_i), Rcpp::Named("j") = Rcpp::wrap(jr),
                            Rcpp::Named("flux") = Rcpp::wrap(res_xr));
};

// Constructeur Residents.
Urban::Residents::Residents(Urban& urb_, int from_) : urb(urb_), from(from_){};

// Constructeurs SubRegion.
Urban::SubRegion::SubRegion(Urban& urb_, std::vector<int> gfrom, std::vector<int> gto)
    : urb(urb_), group_from(gfrom), group_to(gto){};

Urban::SubRegion::SubRegion(Urban& urb_, Rcpp::IntegerVector gfrom, Rcpp::IntegerVector gto)
    : urb(urb_), group_from(Rcpp::as<std::vector<int> >(gfrom)), group_to(Rcpp::as<std::vector<int> >(gto)){};

// Produit le vecteur des colonnes d'emplois disponibles, dans l'ordre des rangs pour les résidents considérées.
std::vector<int> Urban::Residents::map_col_dispo(std::vector<double>& emplois_libres) {
  int debut = urb.p[from];
  int fin = urb.p[from + 1];
  std::vector<int> col_dispo;
  col_dispo.reserve(fin - debut);
  for (auto k = 0; k < fin - debut; ++k) {
    auto pos = urb.jr[debut + k];
    if (emplois_libres[pos] > LIMITE_PRECISION_3) col_dispo.push_back(pos);
  }
  return (col_dispo);
}

// Calcul de l'attractivité des différents sites d'emplois disponible.
// RQ : Le repérage des colonnes disponibles est effectué implicitement.
std::vector<double> Urban::Residents::attractivite(std::vector<double>& emplois_libres, std::shared_ptr<fonction_attraction>& fct) {
  std::vector<int> col_dispo;
  col_dispo = map_col_dispo(emplois_libres);
  int n_sites = col_dispo.size();
  std::vector<double> attract(n_sites);
  for (auto k = 0; k < n_sites; ++k) {
    auto index = col_dispo[k];
    attract[k] = emplois_libres[index] * (*fct)(urb.xr[index]);  // revoir def de la fonction
  }
  return attract;
}

// Calcul de l'attractivité des différents sites d'emplois disponible.
// RQ : Le repérage des colonnes disponibles doit être fourni (pour éviter des recalculs inutiles).
std::vector<double> Urban::Residents::attractivite(std::vector<double>& emplois_libres,
                                                   std::vector<int>& col_dispo,
                                                   std::shared_ptr<fonction_attraction>& fct) {
  int n_sites = col_dispo.size();
  std::vector<double> attract(n_sites);
  for (auto k = 0; k < n_sites; ++k) {
    auto index = col_dispo[k];
    attract[k] = emplois_libres[index] * (*fct)(urb.xr[index]);  // revoir def de la fonction
  }
  return attract;
}

// Calcul la répartition initiale des résidents sur les sites selon leur attractivité.
// i-e la répartition est menée sans égard pour les limitations de places sur chacun des sites
// i-e "ça peut déborder".
// Le mode de calcul est continu.
std::vector<double> Urban::Residents::repartition_nolimit(std::vector<int>& col_dispo, std::vector<double>& attract,
                                                          double actifs_restants) {
  double total = std::accumulate(attract.begin(), attract.end(), 0.0);
  double absorption = -log(urb.fuites[from]) / total;
  int n_sites = attract.size();

  std::vector<double> jobtakers(attract);  // actifs absorbés par sites.
  total = 0.0;
  double anciens, facteur;

  for (auto index = 0; index < n_sites;) {
    auto pos = index + 1;
    while (pos != n_sites && urb.xr[col_dispo[index]] == urb.xr[col_dispo[pos]]) ++pos;
    total = std::accumulate(attract.begin() + index, attract.begin() + pos, 0.0);
    anciens = actifs_restants;
    actifs_restants *= exp(-absorption * total);
    facteur = (anciens - actifs_restants) / total;
    for (auto k = index; k < pos; ++k) {
      jobtakers[k] *= facteur;
    }
    index = pos;
  }
  return (jobtakers);
}

// retourner modifie par mutation le vecteur repartition et renvoie le nombre d'actifs sur le retour.
double retourner(std::vector<double>& repartition, std::vector<double>& places_libres) {
  double retour = 0;
  for (auto k = 0; k < repartition.size(); ++k) {
    if (repartition[k] > places_libres[k]) {
      retour += repartition[k] - places_libres[k];
      repartition[k] = places_libres[k];
    }
  }
  return retour;
}

// Calcule la répartition des résidents sur les sites en tenant compte des limites de places.
// Distribue le nombre d'actifs initiaux (= urb.actifs).
std::vector<double> Urban::Residents::repartition_limited(std::vector<double>& emplois_libres, std::vector<int> col_dispo, std::vector<double> attract) {
  int n_sites = col_dispo.size();
  std::vector<double> repartition(n_sites), places_restantes(n_sites);
  double actifs_restants = urb.actifs[this->from];

  for (auto k = 0; k < n_sites; ++k) places_restantes[k] = emplois_libres[col_dispo[k]];

  double total = std::accumulate(places_restantes.begin(), places_restantes.end(), 0.0);
  if (total < actifs_restants) return places_restantes;

  do {
    // remove-erase idioms
    col_dispo.erase(std::remove_if(col_dispo.begin(), col_dispo.end(),
                                   [&](auto& pos) {
                                     size_t index = &pos - &col_dispo[0];
                                     return places_restantes[index] < LIMITE_PRECISION_3;
                                   }),
                    col_dispo.end());
    attract.erase(std::remove_if(attract.begin(), attract.end(),
                                 [&](auto& pos) {
                                   size_t index = &pos - &attract[0];
                                   return places_restantes[index] < LIMITE_PRECISION_3;
                                 }),
                  attract.end());

    int sites_restants = col_dispo.size();
    std::vector<double> distrib(sites_restants);
    distrib = repartition_nolimit(col_dispo, attract, actifs_restants);
    actifs_restants = retourner(distrib, places_restantes);

    for (auto k = 0, pos = 0; k < sites_restants; ++k, ++pos) {
      while (places_restantes[pos] < LIMITE_PRECISION_3) ++pos;
      repartition[pos] += distrib[k];
      places_restantes[pos] -= distrib[k];
    }
  } while (actifs_restants > LIMITE_PRECISION_3);
  
  return repartition;
}

// regrouper ramasse les flux estimés selon les groupes de SubRegion.
std::vector<float> Urban::SubRegion::regrouper(std::vector<std::vector<float> >& liens) {
  const int N = ng_from();
  const int K = ng_to();
  std::vector<float> flux(N * K);

#pragma omp declare reduction(vsumf : std::vector<float> : std::transform(                    \
        omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<float>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp parallel for reduction(vsumf : flux) collapse(2)
  for (auto i = 0; i < urb.n_from(); ++i) {
    for (auto j = 0; j < urb.n_to(); ++j) {
      flux[group_from[i] * K + group_to[j]] += liens[i][j];
    }
  }
  return flux;
}

// Sortie des flux regroupés.
Rcpp::List Urban::SubRegion::format_sortie(std::vector<std::vector<float> >& liens) {
  const int N = ng_from();
  const int K = ng_to();
  const int taille = N * K;

  std::vector<float> res_xr(taille);
  res_xr = regrouper(liens);

  std::vector<int> res_i(taille);
  std::vector<int> res_j(taille);

  for (auto i = 0; i < N; ++i) {
    for (auto j = 0; j < K; ++j) {
      res_i[i * K + j] = i;
      res_j[i * K + j] = j;
    }
  }
  return Rcpp::List::create(Rcpp::Named("i") = Rcpp::wrap(res_i), Rcpp::Named("j") = Rcpp::wrap(res_j),
                            Rcpp::Named("flux") = Rcpp::wrap(res_xr));
}

// Sortie des flux regroupés et entropie relative avec les flux de référence.
Rcpp::List Urban::SubRegion::format_sortie(std::vector<std::vector<float> >& liens, std::vector<float>& cible) {
  const int N = ng_from();
  const int K = ng_to();
  const int taille = N * K;

  std::vector<float> res_xr(taille);
  res_xr = regrouper(liens);

  std::vector<int> res_i(taille);
  std::vector<int> res_j(taille);

  for (auto i = 0; i < N; ++i) {
    for (auto j = 0; j < K; ++j) {
      res_i[i * K + j] = i;
      res_j[i * K + j] = j;
    }
  }

  float kl = entropie_relative<float>(res_xr, cible);

  return Rcpp::List::create(Rcpp::Named("i") = Rcpp::wrap(res_i), Rcpp::Named("j") = Rcpp::wrap(res_j),
                            Rcpp::Named("flux") = Rcpp::wrap(res_xr), Rcpp::Named("kl") = kl);
}
