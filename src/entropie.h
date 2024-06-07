#ifndef __ENTROPIE__
#define __ENTROPIE__

#include <vector>
#include <numeric>
#include "constants.h"

template <typename T>
T entropie_relative(std::vector<T> estimation, std::vector<T> reference) {
    T tot_estim, tot_ref;
    tot_estim = std::accumulate(estimation.begin(), estimation.end(), 0.0);
    tot_ref = std::accumulate(reference.begin(), reference.end(), 0.0);

    if (tot_estim == 0) Rcpp::stop("Les flux estimés sont tous nuls.");
      
    T kl = 0;
    for (auto k = 0; k < estimation.size(); ++k) {
        estimation[k] /= tot_estim;
        reference[k] /= tot_ref;
        if (reference[k] < PLANCHER_KL) reference[k] = PLANCHER_KL;
        if (estimation[k] != 0) kl += estimation[k] * (log(estimation[k]) - log(reference[k]));
      }

    return kl;
}


#endif //__ENTROPIE__