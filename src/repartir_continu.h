#ifndef __REPARTIR_CONTINU__
#define __REPARTIR_CONTINU__

std::vector<double> one_distrib_continu(const double& entrants, 
                                        const double& fuite,
                                        const std::vector<double>& attractivite,
                                        const std::vector<double>& distances,
                                        const std::vector<double>& placeslibres); 

std::vector<double> repartir_continu(const double actifs, 
                                     const double fuite,
                                     const std::vector<double>& attractivite,
                                     const std::vector<double>& distances,
                                     std::vector<double>& placeslibres);

#endif // __REPARTIR_CONTINU__