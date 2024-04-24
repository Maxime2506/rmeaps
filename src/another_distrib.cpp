#include <Rcpp.h>
using namespace Rcpp;
// Fonction de répartition des actifs entre les sites d'emplois selon l'attractivité du site et la fuite.
// La fonction ne vérifie pas si l'afflux d'actifs excède le nombre de places libres.
// [[Rcpp::export]]
std::vector<double> another_distrib(const double& entrants, 
                                    const double& fuite,
                                    const std::vector<double>& attraction,
                                    const std::vector<double>& xr_dist,
                                    const std::size_t& debut,
                                    const std::vector<double>& placeslibres) {
  
  // Le cas sum(placeslibres) = 0 doit être géré en amont.
  std::size_t k_valid = placeslibres.size();
  std::vector<double> accessibility(k_valid), placesprises(k_valid);
  
  // Calcul de l'accessibilité (corrigée de l'attraction).
  double tot = 0;
  for (std::size_t k = 0; k < k_valid;) {
    auto pos = k + 1L;
    while (xr_dist[debut + k] == xr_dist[debut + pos] && pos < k_valid) ++pos;
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
    while (xr_dist[debut + k] == xr_dist[debut + pos] && pos < k_valid) ++pos;
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