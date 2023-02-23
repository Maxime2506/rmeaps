#ifndef __REPARTIR_ALT__
#define __REPARTIR_ALT__

std::vector<double> repartir_alt(std::vector<double>& placeslibres, 
                                 std::vector<double>& attractivites, 
                                 std::vector<double>& od,
                                 double& fuite, 
                                 double& actifs,
                                 double seuil_newton);

#endif // __REPARTIR_ALT__