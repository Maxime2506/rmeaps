//#include <type_traits>
#include <vector>
#include <cstddef>
#include <iterator>
#include <algorithm>

#include "classes.h"
#include "classes_attraction.h"
#include "constants.h"

Urban::Urban() {};

Urban::Urban(Rcpp::IntegerVector jr_, Rcpp::IntegerVector p_, Rcpp::NumericVector xr_, 
             Rcpp::NumericVector actifs_, Rcpp::NumericVector emplois_, Rcpp::NumericVector fuites_):
            jr(Rcpp::as< std::vector<unsigned int> >(jr_)), 
            p(Rcpp::as< std::vector<unsigned long> >(p_)), 
            xr(Rcpp::as< std::vector<double> >(xr_)), 
            actifs(Rcpp::as< std::vector<double> >(actifs_)), 
            emplois(Rcpp::as< std::vector<double> >(emplois_)), 
            fuites(Rcpp::as< std::vector<double> >(fuites_)) {
  actifs_libres.resize(actifs_.size());
  emplois_libres.resize(emplois_.size());
            std::copy(actifs.begin(), actifs.end(), actifs_libres.begin()); 
            std::copy(emplois.begin(), emplois.end(), emplois_libres.begin());
};

Urban::Urban(const std::vector<unsigned int> jr_, const std::vector<unsigned long> p_, const std::vector<double> xr_, 
      const std::vector<double> actifs_, const std::vector<double> emplois_, const std::vector<double> fuites_):
    jr(jr_), p(p_), xr(xr_), actifs(actifs_), emplois(emplois_), fuites(fuites_) {
  actifs_libres.resize(actifs_.size());
  emplois_libres.resize(emplois_.size());
    std::copy(actifs.begin(), actifs.end(), actifs_libres.begin()); 
    std::copy(emplois.begin(), emplois.end(), emplois_libres.begin());
};

Urban::Residents::Residents(Urban& urbs_, unsigned int from_): 
    urbs(urbs_), from(from_), col_dispo(map_col_dispo(urbs_, from_)) {};

std::vector<unsigned int> map_col_dispo(Urban& urbs_, unsigned int from_) {
    unsigned long debut = urbs_.p[from_];
    unsigned long fin = urbs_.p[from_ + 1];
    std::vector<unsigned int> col_dispo;
    col_dispo.reserve(fin - debut);
    for(auto k = 0; k < fin - debut; ++k) {
        auto pos = urbs_.jr[debut + k];
        if (urbs_.emplois_libres[pos] > LIMITE_PRECISION_3) col_dispo.push_back(pos);
    }
    return(col_dispo);
}

std::vector<double> attractivite(Urban::Residents& res, 
                                 std::shared_ptr<fonction_attraction>& fct) {
    // Calcul de l'attractivité.
    unsigned int n_sites = res.col_dispo.size();
    std::vector<double> attract(n_sites);
    for (auto k = 0; k < n_sites; ++k) {
        auto index = res.col_dispo[k];
        attract[k] = res.urbs.emplois_libres[index] * (*fct)(res.urbs.xr[index]); // revoir def de la fonction
    }
    return attract;
}
std::vector<double> attract_cumul(Urban::Residents& res, 
                                  std::vector<double>& attract) {
    // Calcul du cumul avec la même valeur pour les équidistants.
    std::vector<double> cumul(attract.size());
    double tot = 0;
    for (auto k = res.col_dispo.begin(); k != res.col_dispo.end();) {
        auto pos = std::next(k);
        while (pos != res.col_dispo.end() && res.urbs.xr[*k] == res.urbs.xr[*pos]) ++pos;
        for (auto ego = k; ego != pos; ++ego) {
            auto index = std::distance(res.col_dispo.begin(), ego);
            tot += attract[index];
        }
        for (auto index = std::distance(res.col_dispo.begin(), k); index < std::distance(res.col_dispo.begin(), pos); ++index) {
            cumul[index] = tot;
        }
        k = pos;
    }
    return cumul;
}

std::vector<double> calc_repartition(Urban::Residents& res, 
                                     std::vector<double>& attract,
                                     std::vector<double>& cumul) {
    // Calcul de l'absorption sur la ligne from considérée.
    double absorption = -log(res.urbs.fuites[res.from]) / cumul.back();
    unsigned int n_sites = cumul.size();

    // Calcul des actifs absorbés par sites.
    std::vector<double> jobtakers(n_sites + 1, res.urbs.actifs_libres[res.from]);
    for(auto k = 0; k < n_sites; ++k) {
        jobtakers[k + 1] *= exp(-absorption * cumul[k]); // ceux qui dépassent le site k+1.
    }
    for(auto k = 0; k < n_sites; ++k) {
        jobtakers[k] -= jobtakers[k + 1];
    }
    // Répartition des jobtakers en traitant les cas à distances égales.
    auto debut = res.col_dispo.begin();
    for (auto k = debut; k != res.col_dispo.end();) {
        double tot_attract = 0, tot_jobtakers = 0;
        auto pos = std::next(k);
        while (pos != res.col_dispo.end() && res.urbs.xr[*k] == res.urbs.xr[*pos]) ++pos;
        for (auto it = k; it != pos; ++it) {
            auto index = std::distance(debut, it);
            tot_attract += attract[index];
            tot_jobtakers += jobtakers[index];
        }
        if (tot_attract > 0) {
            for (auto it = k; it != pos; ++it) {
                auto index = std::distance(debut, it);
                jobtakers[index] = attract[index] / tot_attract * tot_jobtakers;
            }
        }
        k = pos;
    }
    jobtakers.resize(n_sites);
    return jobtakers;
}
