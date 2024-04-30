#ifndef __MEAPS_CORE__
#define __MEAPS_CORE__


std::vector< std::vector<double> > meaps_core(const std::vector<int> jr_dist, 
                                              const std::vector<int> p_dist, 
                                              const std::vector<double> xr_dist, 
                                              std::vector<double> emplois,
                                              const std::vector<double> actifs, 
                                              std::vector<double> fuite, 
                                              const std::vector<double> parametres,
                                              const std::vector<double> xr_odds,
                                              const std::string attraction,
                                              const int nthreads, 
                                              const bool verbose);


#endif // __MEAPS_CORE__