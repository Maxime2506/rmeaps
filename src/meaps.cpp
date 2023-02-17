#include <Rcpp.h>
#include "repartir_actifs.h"
using namespace Rcpp;
using namespace std;

//' La fonction meaps sur un shuf. 
//' @param rkdist La matrice des rangs dans lequel les colonnes j sont passées en revue pour chacune des lignes i.
//' @param emplois Le vecteur des emplois disponibles sur chacun des sites j (= marge des colonnes).
//' @param actifs Le vecteur des actifs partant de chacune des lignes visées par shuf. 
//' @param modds La matrice des odds modifiant la chance d'absorption de chacun des sites j pour des résidents en i.
//' @param f Le vecteur de la probabilité de fuite des actifs hors de la zone d'étude.
//' @param shuf Le vecteur de priorité des actifs pour choisir leur site d'arrivée. 
//'        Il est possible de segmenter les départs d'une ligne i en répétant cette ligne à plusieurs endroits du shuf.
//'        Dans ce cas, le nombre d'actifs sera répartie également entre les différents départs depuis cette ligne.
//' @param fuite_min Seuil minimal pour la fuite d'un actif. Doit être supérieur à 0. Défault = 1e-3.
//' @param seuil_newton Seuil de convergence pour la méthde de Newton du calcul des probas d'absorption.
//' 
//' @return renvoie une matrice avec les estimations du nombre de trajets de i vers j.
// [[Rcpp::export]]
NumericMatrix meaps_oneshuf(
    const IntegerMatrix rkdist, 
    NumericVector emplois,
    const NumericVector actifs,
    const NumericMatrix modds,
    const NumericVector f, 
    const IntegerVector shuf,
    bool normalisation = false,
    double fuite_min = 1e-3,
    double seuil_newton = 1e-6)
{
  const int N = rkdist.nrow(),
      K = rkdist.ncol(),
      Ns = shuf.size();
  int k_valid;
  NumericMatrix liaisons(N,K);
  NumericVector odds(K);
  IntegerVector rki(K), ishuf(Ns);
  LogicalVector nna_rki(K);
  
  double tot_emp = sum(emplois);
  double tot_act_inzone = sum(actifs * (1 - f));
  // Attention : calage des emplois sur le nombre d'actifs.
  if (normalisation) {
      emplois = emplois * tot_act_inzone / tot_emp;
  }
  vector<double> emp(K); 
  emp = as<vector<double>>(emplois);
  ishuf = shuf - 1L;
  
  // Le vecteur shuf peut être plus long que le nombre de départs d'actifs s'il fait repasser plusieurs fois
  // la même ligne d'actifs. Dans ce cas, on compte la fréquence de passage de chaque ligne et l'on divise le
  // poids de la ligne par cette fréquence.
  std::vector<int> freq_actifs(N, 0L);
  for (auto i: ishuf) {
    freq_actifs[i]++;
  }
  
  for (int i: ishuf) {
    
    // On vérifie si on est pas trop long
    if(i%1000 == 1) {Rcpp::checkUserInterrupt();}

    rki = rkdist(i, _);
    nna_rki = !is_na(rki);
    k_valid = sum(!is_na(rki));
    IntegerVector arrangement(k_valid);
    
    for(int j = 0; j < K; ++j) {
      if(nna_rki[j]==TRUE)
        arrangement[rki[j] - 1L] = j;
    }
    
    vector<double> dispo (k_valid), od(k_valid), repartition(k_valid);
    for (int j = 0; j < k_valid; ++j) {
      dispo[j] = emp[arrangement[j]];
      od[j] = modds(i, arrangement[j]);
    }
    
    // Choix d'une limite basse pour la fuite.
    double fuite = max(fuite_min, f[i]);
    double actifsinzone = (1 - fuite) * actifs[i] / freq_actifs[i];
    
    repartition = repartir_actifs(dispo, od, fuite, actifsinzone, seuil_newton);
    
    // Inscription des résultats locaux dans la perspective globale.
    for(int j = 0; j < k_valid ; ++j) {
      emp[arrangement[j]] -= repartition[j];
      liaisons(i, arrangement[j]) += repartition[j];
    }
  }
  return liaisons;
}

