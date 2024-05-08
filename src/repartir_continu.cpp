
#include <algorithm>
#include <iterator>
#include <vector>
#include <cmath>

std::vector<double> one_distrib_continu(const double entrants, 
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
std::vector<double> repartir_continu(const double actifs, 
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
    nouvellesprises = one_distrib_continu(actifs_non_etablis, fuite, attractivite, distances, placesrestantes);
    
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