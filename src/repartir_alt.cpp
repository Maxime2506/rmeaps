#include <Rcpp.h>
#include "fonctions_newton.h"

using namespace Rcpp;

std::vector<double> repartir_alt(std::vector<double>& placeslibres, 
                                 std::vector<double>& attractivites,
                                 std::vector<double>& od,
                                 double& fuite, 
                                 double& actifs,
                                 double seuil_newton = 1e-6) {
  
  int k_valid = placeslibres.size();
  
  // Calcul de la place disponible total sur la ligne.
  double tot = 0;
  for (auto& d: placeslibres) {
    tot += d;
  }
  // rq : lorsque tot < actifs restants, tous les emplois dispo sont pris.
  if (tot <= actifs * (1 - fuite)) {
    return placeslibres; // c'est ici qu'il peut y avoir une fuite supplémentaire.
  }
  
  std::vector<double> repartition(k_valid, 0.0);
  // if (k_valid == 1) {
  //   repartition[0] = std::min(dispo[0], actifs);
  //   return repartition;
  // }
  
  double p_ref, c_ref, new_cref, eps;
  p_ref = 1  - pow(fuite, 1 / tot);
  // Calcul par la méthode de Newton de la chance d'absorption de référence compatible avec la fuite.
  c_ref = p_ref / (1 - p_ref);// Chance d'absorption de référence. Calcul initial non calé sur la fuite.
  
  int compteur = 0;
  do {
    new_cref = c_ref - (sumlog_passage(c_ref, placeslibres, attractivites, od) + log(fuite))/ d_sumlog_passage(c_ref, placeslibres, attractivites, od); 
    eps = std::abs(new_cref - c_ref);
    c_ref = new_cref;
    compteur++;
    if (compteur > 100) {
      Rcout << "Le calcul par la méthode de Newton de la chance d\'absorption n\'a pas convergé.\n";
      eps = 0; // L'estimation est poursuivie malgré tout...
    }
  } while (eps > seuil_newton);
   
  // Calcul des actifs arrivants sur chacun des sites d'emplois, en tenant compte de la capacité en regard de la demande (débordement).
  double actifsrestants = actifs, debordement = 0.0, demande;
  for (int j = 0; j < k_valid; ++j) {
    if (placeslibres[j] > 0.0) {
      demande = actifsrestants * (1 - pow(1 + od[j] * c_ref, -attractivites[j]));
      if (demande > placeslibres[j]) {
        repartition[j] = placeslibres[j];
        debordement += demande - placeslibres[j];
      } else {
        repartition[j] = demande;
      }
      actifsrestants -= demande;
    }
  }
  
  // Replacement des débordements sur les sites non remplis à proportion de la répartition calculée ci-dessus.``
  // RQ : debordement peut atteindre 0 car on a évacué dès le départ le cas où il n'y a pas assez d'emplois pour le nombre d'actifs.
  double repartition_tot, new_debord, demande_sup; 
  while (debordement > 0) {
    repartition_tot = 0.0;
    new_debord = 0.0;

    for (int j = 0; j < k_valid; ++j) {
      if (repartition[j] < placeslibres[j]) { repartition_tot += repartition[j]; }
    }
    
    for (int j = 0; j < k_valid; ++j)
      { if (placeslibres[j] > 0.0) {
        if (repartition[j] < placeslibres[j]) {
          demande_sup = debordement * repartition[j]  / repartition_tot;
          if (demande_sup + repartition[j] > placeslibres[j]) {
            new_debord += repartition[j] + demande_sup - placeslibres[j];
            repartition[j] = placeslibres[j];
          } else {
            repartition[j] += demande_sup; 
          }
        }
      }
      debordement = new_debord;
      }
  }
  return repartition;
} 