#ifndef __REPARTIR_ACTIFS__
#define __REPARTIR_ACTIFS__

std::vector<double> repartir_actifs(std::vector<double>& dispo, 
                                    std::vector<double>& od,
                                    double& fuite, 
                                    double& actifs,
                                    double seuil_newton);

#endif // __REPARTIR_ACTIFS__