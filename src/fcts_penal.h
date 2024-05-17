#ifndef __FCTS_PENAL__
#define __FCTS_PENAL__

 double marche(double x, const double rayon, const double plafond);
 
 double marche_liss(double x, const double rayon, const double plafond);

 double double_marche_liss(double x, const double r1, const double r2, const double o1, const double o2);

 double decay(double x, const double delta, const double plancher);
 
 double logistique(double x, const double rayon, const double amplitude, const double plancher);

#endif // __FCTS_PENAL__