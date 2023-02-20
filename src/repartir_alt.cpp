#include <Rcpp.h>
#include "deborder.h"
#include "distribuer.h"
#include "utils_newton_methods.h"

using namespace Rcpp;

std::vector<double> repartir_alt(std::vector<double>& dispo, 
                                 std::vector<double>& od,
                                 double& fuite, 
                                 double& actifs,
                                 double seuil_newton = 1e-6) {
  
  int k_valid = dispo.size();
  std::vector<double> repartition(k_valid);
  
  // Calcul de la place disponible total sur la ligne.
  double tot = 0;
  for (auto& d: dispo) {
    tot += d;
  }
  if (tot < 1e-2 || k_valid == 1) {
    // Seuil où l'emploi disponible est suffisamment petit pour s'épargner des calculs inutiles et fragiles.
    repartition = deborder(dispo, actifs);
    return repartition ;
  }
  
  std::vector<double> proportions(k_valid);
  double p_ref, c_ref;
  p_ref = 1  - pow(fuite, 1 / tot);
  
  if (p_ref < 0.1) { // seuil de référence en deçà duquel on effectue un calcul approché des log.
    
    // Calcul approché de c_ref en approchant log[(1 + c_abs)^dispo] ~ dispo * od * c_ref
    double temp = 0.0;
    for (int j = 0; j < k_valid; ++j) {
      temp += dispo[j] * od[j];
    }
    c_ref = - log(fuite) / temp;
    
  } else {
    
    // Calcul par la méthode de Newton de la chance d'absorption de référence compatible avec la fuite.
    double new_cref, eps;
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
  }
  
  // Calcul des proba d'arrivées sur chacun des sites (dépend du chemin, != p_abs).
  double logpass = 0.0, logfuit;
  for (int j = 0; j < k_valid; ++j) {
    logfuit = - dispo[j] * log(1 + od[j] * c_ref); // proba conditionnelle en log de fuir j une fois arrivée jusqu'à j.
    proportions[j] = exp(logpass) * (1 - exp(logfuit));
    logpass += logfuit; // proba en log d'arriver jusqu'à j+1, calculé à partir de la proba d'arriver jusqu'en j.
  }
  
  repartition = distribuer(dispo, proportions, actifs);
  
  return repartition;
} 