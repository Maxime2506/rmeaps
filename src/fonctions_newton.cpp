#include <Rcpp.h>
using namespace std;

//' Fonction de calcul du logarithme des probabilités de fuir.
//' @param x Cette valeur renvoie à c_ref = la chance d'être absorbé par un site sous l'hypothèse d'homogénéité.
//' @param placeslibres La capacité de chacun des sites.
//' @param attractivites L'attractivité de chacun des sites (puissance appliqué à la proba de passer au site suivant).
//' @param odds Les odds modifiant l'absorption relative de chacun des sites.
//' 
//' @return La fonction renvoie le log de la fuite.
double sumlog_passage(
    double &x,
    vector<double>& placeslibres,
    vector<double>& attractivites,
    vector<double>& odds) {
  
  int n = placeslibres.size();
  double res = 0.;
  
  for(int i = 0; i < n; ++i) {
    if (placeslibres[i] > 0) {
      res += attractivites[i] * log(1. + odds[i] * x);
      }
    }
  
  return (res);
}

//' Fonction de calcul de la dérivée de la fonction log_fuite.
//' @param x Cette valeur renvoie à c_ref = la chance d'être absorbé par un site sous l'hypothèse d'homogénéité.
//' @param placeslibres La capacité de chacun des sites.
//' @param attractivites L'attractivité de chacun des sites (puissance appliqué à la proba de passer au site suivant).//' @param odds Les odds modifiant l'absorption relative de chacun des sites.
//' 
//' @return La fonction renvoie la dérivée du logarithme de la fonction de fuite.
double d_sumlog_passage(
    double &x,
    vector<double>& placeslibres,
    vector<double>& attractivites,
    vector<double>& odds) {
  
  int n = placeslibres.size();
  double res = 0.;
  
  for(int i = 0; i < n; ++i) {
    if (placeslibres[i] > 0) {
      res += attractivites[i] * odds[i] / (1. + odds[i] * x) ;
      }
    }
  return (res);
}
