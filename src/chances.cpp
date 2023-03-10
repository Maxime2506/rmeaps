#include <Rcpp.h>
#include "fonctions_newton.h"
using namespace Rcpp;
using namespace std;

//' La fonction renvoie les chances d'absorption issue de meaps pour la méthode continue (default). 
//' @param rkdist La matrice des rangs dans lequel les colonnes j sont passées en revue pour chacune des lignes i.
//' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes).
//' @param modds La matrice des odds modifiant la chance d'absorption de chacun des sites j pour des résidents en i.
//' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude.
//' 
//' @return renvoie le vecteur des chances d'absorption.
// [[Rcpp::export]]
NumericMatrix chances_absorption(
    const IntegerMatrix rkdist, 
    const NumericVector emplois,
    const NumericMatrix modds,
    const NumericVector f)
{
  const int N = rkdist.nrow(),
            K = rkdist.ncol();
  int k_valid;
  double tot, p_ref, c_ref, new_cref, eps;
  
  NumericVector odds(K);
  IntegerVector rki(K);
  LogicalVector nna_rki(K);
  NumericMatrix chances(N, K);
  
  tot = 0;
  for (auto& k: emplois) {
    tot += k;
  }
  
  for (int i = 0; i < N; ++i) {
    
    rki = rkdist(i, _);
    nna_rki = !is_na(rki);
    k_valid = sum(!is_na(rki));
    IntegerVector arrangement(k_valid);
    
    for(int j = 0; j < K; ++j) {
      if(nna_rki[j]==TRUE)
        arrangement[ rki[j] - 1L ] = j;
    } 
    
    vector<double> dispo (k_valid);
    for (int j = 0; j < k_valid; ++j) {
      dispo[j] = emplois[ arrangement[j] ]; 
    }
    
    // Choix d'une limite basse pour la fuite.
    double fuite = max(1e-3, f[i]);
    
    odds = modds(i, _);
    p_ref = 1  - pow(fuite, 1 / tot);
    c_ref = p_ref / (1 - p_ref);// Chance d'absorption de référence. Calcul initial non calé sur la fuite.
    
    vector<double> od = as<vector<double>>(odds[arrangement]);
    // Calcul par la méthode de Newton de la chance d'absorption de référence compatible avec la fuite.
    do {
      new_cref = c_ref - (sumlog_passage(c_ref, dispo, dispo, od) + log(fuite))/ d_sumlog_passage(c_ref, dispo, dispo, od); 
      eps = abs(new_cref - c_ref);
      c_ref = std::max(0., new_cref);
    } while (eps > 1e-6);
    
    vector<double> c_abs(k_valid);
    for (int j = 0; j < k_valid; ++j) {
      chances(i, arrangement[j]) = od[j] * c_ref; // chance d'absorption.
    }
  }

return chances;
}

