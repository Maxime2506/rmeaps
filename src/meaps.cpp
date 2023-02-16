#include <Rcpp.h>
#include "deborder.h"
#include "distribuer.h"
#include "utils_newton_methods.h"
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
//' 
//' @return renvoie une matrice avec les estimations du nombre de trajets de i vers j.
// [[Rcpp::export]]
NumericMatrix meaps_oneshuf(
    const IntegerMatrix rkdist, 
    const NumericVector emplois,
    const NumericVector actifs,
    const NumericMatrix modds,
    const NumericVector f, 
    IntegerVector shuf)
{
  int N = rkdist.nrow(),
      K = rkdist.ncol(),
      Ns = shuf.size(),
      k_valid;
  double tot, p_ref, c_ref, new_cref, eps;
  NumericMatrix liaisons(N,K);
  NumericVector odds(K);
  IntegerVector rki(K), ishuf(Ns);
  LogicalVector nna_rki(K);
  
  vector<double> emp(K); 
  // Attention : calage des emplois sur le nombre d'actifs.
  double tot_emp = sum(emplois);
  double tot_act_inzone = sum(actifs * (1 - f));
  for(int j = 0; j < K; ++j) {
    emp[j] = tot_act_inzone * emplois[j] / tot_emp;
  }
  
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
    // Si aucun actif au départ, ou valeur négative, on ne calcule rien.

    rki = rkdist(i, _);
    nna_rki = !is_na(rki);
    k_valid = sum(!is_na(rki));
    IntegerVector arrangement(k_valid);
    
    for(int j = 0; j < K; ++j) {
      if(nna_rki[j]==TRUE)
        arrangement[rki[j]-1L] = j;
    }
    
    vector<double> dispo (k_valid), repartition(k_valid);
    for (int j = 0; j < k_valid; ++j) {
      dispo[j] = emp[arrangement[j]]; 
    }

    tot = 0;
    for (auto& d: dispo) {
      tot += d;
      }
    // Choix d'une limite basse pour la fuite.
    double fuite = max(1e-3, f[i]);
    
    // Seuil où l'emploi disponible est suffisamment petit pour s'épargner des calculs inutiles et fragiles.
    if (tot < 1e-2 || k_valid == 1) {
      
      repartition = deborder(dispo, (1 - fuite) * actifs[i] / freq_actifs[i]);
      
    } else {
      odds = modds(i, _);
      p_ref = 1  - pow(fuite, 1 / tot);
      c_ref = p_ref / (1 - p_ref);// Chance d'absorption de référence. Calcul initial non calé sur la fuite.
      
      vector<double> od = as<vector<double>>(odds[arrangement]),
                     proportions(k_valid);
      // Calcul par la méthode de Newton de la chance d'absorption de référence compatible avec la fuite.
      do {
        new_cref = c_ref - (log_fuite(c_ref, dispo, od) + log(fuite))/ d_logfuite(c_ref, dispo, od); 
        eps = abs(new_cref - c_ref);
        c_ref = new_cref;
      } while (eps > 1e-6);
      
      vector<double> c_abs(k_valid);
      for (int j = 0; j < k_valid; ++j) {
        c_abs[j] = od[j] * c_ref;
        }
      
      double logpass = 0.0, 
             logfuit;
      
      for (int j = 0; j < k_valid; ++j) {
        logfuit = - dispo[j] * log(1 + c_abs[j]); // proba conditionnelle en log de fuir j une fois arrivée jusqu'à j.
        proportions[j] = exp(logpass) * (1 - exp(logfuit));
        logpass += logfuit; // proba en log d'arriver jusqu'à j+1, calculé à partir de la proba d'arriver jusqu'en j.
      }
      repartition = distribuer(dispo, proportions, (1 - fuite) * actifs[i] / freq_actifs[i]);
   }
    
    // Inscription des résultats locaux dans la perspective globale.
    for(int j = 0; j < k_valid ; ++j) {
      emp[arrangement[j]] -= repartition[j];
      liaisons(i, arrangement[j]) += repartition[j];
    }
  }
  return liaisons;
}

