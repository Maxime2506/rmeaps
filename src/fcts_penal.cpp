// Fonctions de pénalités utilisées pour modifier l'attractivité d'un site dans MEAPS.
// Remarque : il suffit de définir ces fonctions à un facteur près en raison d'un renormalisation implicite dans MEAPS.
#include <cmath>

//' Fonction de pénalité "marche" : vaut 1 sur un rayon fixé, et decru au-delà.
//' @param x distance.
//' @param rayon distance de la marche.
//' @param plancher point bas après la marche.
//' 
//' @return un facteur d'attraction
 double marche(double x, const double rayon, const double plafond) {
   if (x > rayon) return 1;
   return plafond;
 }
 
//' Fonction de pénalité "marche_liss" : vaut 1 sur un rayon fixé, et décroit au-delà,
//' en lissant le passge pour les distances fractionaires
//' @param x distance.
//' @param rayon distance de la marche.
//' @param plancher point bas après la marche.
//' 
//' @return un facteur d'attraction
 double marche_liss(double x, const double rayon, const double plafond) {
   if (x > ceil(rayon)) return 1;
   if (x <= floor(rayon)) return plafond;
   return plafond + (1 - plafond) *  (ceil(rayon) - rayon);
 }
 
//' Fonction de pénalité "double_marche_liss" : vaut 1 au delà de r1+r2+1, et prend deux valeurs (o1, o2),
//' lorsque x<=r1 et r1<x<=r2+r1,
//' en lissant le passage pour les distances fractionaires
//' @param x distance.
//' @param rayon distance de la marche.
//' @param plancher point bas après la marche.
//' 
//' @return un facteur d'attraction
 double marche_liss(double x, const double r1, const double r2,
                    const double o1, const double o2) {
   if (x > ceil(r1+r2)) return 1;
   if (x <= floor(r1)) return o1;
   if (x < ceil(r1)) return o1 + (o2 - o1) *  (ceil(r1) - r1);
   if (x <= floor(r1+r2)) return o2;
   return o2 + (1 - o2) *  (ceil(r1+r2) - r1 - r2);
 }

//' Fonction de pénalité "decay" : odd = 1/d^delta + plancher,
//' @param x distance.
//' @param delta, exposant de la distance.
//' @param plancher valeur pour les distances infinies.
//' 
//' @return un facteur d'attraction
 double decay(double x, const double delta, const double plancher) {
   return exp(-delta * log(x)) + plancher;
 }
 
// Fonction de type logistique x -> 1 + amplitude * { exp(-(x-rayon)) - 1} / { exp(-(x-rayon)) + 1 }
//' @param x distance.
//' @param rayon distance de la bascule.
//' @param amplitude raideur de la bascule.
//' @param plancher point bas après la marche.
//' 
//' @return un facteur d'attraction
 double logistique(double x, const double rayon, const double amplitude, const double plancher) {
   double ex = exp( (rayon-x)/amplitude );
   return plancher + ex / (ex + 1);
 }
