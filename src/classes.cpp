//#include <type_traits>
#include <vector>

#include "classes.h"
#include "classes_attraction.h"
#include "constants.h"

Urban::Urban() { 
    std::vector<unsigned int> jr_(1); 
    std::vector<std::size_t> p_(1);
    std::vector<double> xr_(1);
    std::vector<double>& actifs_(1); 
    std::vector<double>& emplois_(1); 
    std::vector<double>& fuites_(1);
    std::vector<double> actifs_libres_(1);
    std::vector<double> emplois_libres_(1);
    jr = jr_; p = p_; xr = xr_; actifs = actifs_; emplois = emplois_; fuite = fuite_; 
    actifs_libres = actifs_libres_; emplois_libres = emplois_libres_;
}

Urban::Urban(std::vector<unsigned int>& jr_, std::vector<std::size_t>& p_, std::vector<double>& xr_, 
    std::vector<double>& actifs_, std::vector<double>& emplois_, std::vector<double>& fuites_,
    std::vector<double> actifs_libres_, std::vector<double> emplois_libres_):
        jr(jr_), p(p_), xr(xr_), actifs(actifs_), emplois(emplois_), fuites(fuites_),
        actifs_libres(actifs_libres_), emplois_libres(emplois_libres_) {};

Urban::Residents::Residents(Urban& urbs_, unsigned int from_): 
    urbs(urbs_), from(from_), col_dispo(map_col_dispo(urbs_, from_)) {};

std::vector<unsigned int> map_col_dispo(Urban& urbs_, unsigned int from_) {
    std::size_t debut = urbs_.p[from_];
    std::size_t fin = urbs_.p[from_ + 1];
    col_dispo.reserve(fin - debut);
    for(auto k = 0; k < fin - debut; ++k) {
        auto pos = urbs_.jr[debut + k];
        if (urbs_.emplois_libres[pos] > LIMITE_PRECISION_3) col_dispo.push_back(pos);
    }
    return(col_dispo);
}

std::vector<double> attractivite(Urban::Residents& res(urbs_, from_), 
                                 std::unique_ptr<fonction_attraction> fct) {
    // Calcul de l'attractivité.
    std::vector<double> attract(res.col_dispo.size());
    for (auto k: res.col_dispo) {
        auto index = std::distance(res.col_dispo.begin(), k);
        attract[index] = urbs_.emplois_libres[*k] * fct(urbs_.xr[*k]); // revoir def de la fonction
    }
    return attract;
}
std::vector<double> attract_cumul(Urban::Residents& res(urbs_, from_), 
                                  std::vector<double>& attract) {
    // Calcul du cumul avec la même valeur pour les équidistants.
    std::vector<double> cumul(attract.size());
    double tot = 0;
    for (auto k = res.col_dispo.begin(); k != res.col_dispo.end();) {
        auto pos = std::next(k);
        while (pos != res.col_dispo.end() && urbs_.xr[*k] == urbs_.xr[*pos]) ++pos;
        for (auto ego = k; ego != pos; ++ego) {
            auto index = std::distance(res.col_dispo.begin(), ego);
            tot += attract[index];
        }
        for (auto index = std::distance(res.col_dispo.begin(), k); index < std::distance(res.col_dispo.begin(), ego); ++index) {
            cumul[index] = tot;
        }
        k = pos;
    }
    return cumul;
}

std::vector<double> calc_repartition(Urban::Residents& res(urbs_, from_), 
                                     std::vector<double>& attract,
                                     std::vector<double>& cumul) {
    // Calcul de l'absorption sur la ligne from considérée.
    double absorption = -log(urbs_.fuite[from_]) / cumul.back();
    unsigned int n_sites = cumul.size();

    // Calcul des actifs absorbés par sites.
    std::vector<double> jobtakers(n_sites + 1L, urbs_.actifs_libres[from_]);
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
        while (pos != res.col_dispo.end() && urbs_.xr[*k] == urbs_.xr[*pos]) ++pos;
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
