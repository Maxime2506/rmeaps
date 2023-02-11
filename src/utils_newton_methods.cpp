#include <Rcpp.h>
using namespace std;

//' Fonction de calcul du logarithme des probabilités de fuir.
//' @param x Cette valeur renvoie à c_ref = la chance d'être absorbé par un site sous l'hypothèse d'homogénéité.
//' @param dispo La taille de chacun des sites.
//' @param odds Les odds modifiant l'absorption relative de chacun des sites.
//' 
//' @return La fonction renvoie le log de la fuite.
double log_fuite(
    double &x,
    vector<double>& dispo,
    vector<double>& odds) {
  
  int n = dispo.size();
  double res = 0.;
  
  for(int i = 0; i < n; ++i) res += dispo[i] * log(1. + odds[i] * x);
  
  return (res);
}

//' Fonction de calcul de la dérivée de la fonction log_fuite.
//' @param x Cette valeur renvoie à c_ref = la chance d'être absorbé par un site sous l'hypothèse d'homogénéité.
//' @param dispo La taille de chacun des sites.
//' @param odds Les odds modifiant l'absorption relative de chacun des sites.
//' 
//' @return La fonction renvoie la dérivée du logarithme de la fonction de fuite.
double d_logfuite(
    double &x,
    vector<double>& dispo,
    vector<double>& odds) {
  
  int n = dispo.size();
  double res = 0.;
  
  for(int i = 0; i < n; ++i) res += dispo[i] * odds[i] / (1. + odds[i] * x) ;
  
  return (res);
}
