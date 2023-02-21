#include <Rcpp.h>
#include "utils_newton_methods.h"

using namespace Rcpp;

std::vector<double> repartir_alt(std::vector<double>& dispo, 
                                 std::vector<double>& od,
                                 double& fuite, 
                                 double& actifs,
                                 double seuil_newton = 1e-6) {
  
  int k_valid = dispo.size();
  
  // Calcul de la place disponible total sur la ligne.
  double tot = 0;
  for (auto& d: dispo) {
    tot += d;
  }
  // rq : lorsque tot < actifs, tous les emplois dispo sont pris.
  if (tot <= actifs) {
    return dispo;
  }
  
  std::vector<double> repartition(k_valid);
  if (k_valid == 1) {
    repartition[0] = std::min(dispo[0], actifs);
    return repartition;
  }
  
  double p_ref, c_ref, new_cref, eps;
  p_ref = 1  - pow(fuite, 1 / tot);
  // Calcul par la méthode de Newton de la chance d'absorption de référence compatible avec la fuite.
  c_ref = p_ref / (1 - p_ref);// Chance d'absorption de référence. Calcul initial non calé sur la fuite.
  
  int compteur = 0;
  do {
    new_cref = c_ref - (log_fuite(c_ref, dispo, od) + log(fuite))/ d_logfuite(c_ref, dispo, od); 
    eps = std::abs(new_cref - c_ref);
    c_ref = new_cref;
    compteur++;
    if (compteur > 200) {
      Rcout << "Le calcul par la méthode de Newton de la chance d\'absorption n\'a pas convergé.\n";
      eps = 0; // L'estimation est poursuivie malgré tout...
    }
  } while (eps > seuil_newton);
  
  // Calcul des proba d'arrivées sur chacun des sites (dépend du chemin, != p_abs).
  double actifsrestants = actifs;
  for (int j = 0; j < k_valid; ++j) {
    repartition[j] = actifsrestants * (1 - pow(1 + od[j] * c_ref, -dispo[j]));
    actifsrestants -= repartition[j];
    }
  
  return repartition;
} 