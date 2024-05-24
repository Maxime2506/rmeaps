#include <vector>
#include <iterator>
#include <algorithm>

#include "classes.h"
#include "classes_attraction.h"
#include "constants.h"
#include "entropie.h"

Urban::Urban() {};

Urban::Urban(Rcpp::IntegerVector jr_, Rcpp::IntegerVector p_, Rcpp::NumericVector xr_, 
             Rcpp::NumericVector actifs_, Rcpp::NumericVector emplois_, Rcpp::NumericVector fuites_):
            jr(Rcpp::as< std::vector<int> >(jr_)), 
            p(Rcpp::as< std::vector<int> >(p_)), 
            xr(Rcpp::as< std::vector<double> >(xr_)), 
            actifs(Rcpp::as< std::vector<double> >(actifs_)), 
            emplois(Rcpp::as< std::vector<double> >(emplois_)), 
            fuites(Rcpp::as< std::vector<double> >(fuites_)) {
  actifs_libres.resize(actifs_.size());
  emplois_libres.resize(emplois_.size());
            std::copy(actifs.begin(), actifs.end(), actifs_libres.begin()); 
            std::copy(emplois.begin(), emplois.end(), emplois_libres.begin());
};

Urban::Urban(const std::vector<int> jr_, const std::vector<int> p_, const std::vector<double> xr_, 
      const std::vector<double> actifs_, const std::vector<double> emplois_, const std::vector<double> fuites_):
    jr(jr_), p(p_), xr(xr_), actifs(actifs_), emplois(emplois_), fuites(fuites_) {
  actifs_libres.resize(actifs_.size());
  emplois_libres.resize(emplois_.size());
    std::copy(actifs.begin(), actifs.end(), actifs_libres.begin()); 
    std::copy(emplois.begin(), emplois.end(), emplois_libres.begin());
};

Rcpp::List Urban::procurer(std::vector< std::vector<float> >& liens) {
  
    std::vector<float> res_xr(xr.size());
    std::vector<int> res_i(xr.size());

#pragma omp parallel for
    for (auto i = 0; i < n_from(); ++i) {
        for (auto k = p[i]; k < p[i + 1]; ++k) {
            res_xr[k] = liens[i][ jr[k] ];
            res_i[k] = i;
        }
    }
    return Rcpp::List::create( Rcpp::Named("i") = Rcpp::wrap(res_i), 
                               Rcpp::Named("j") = Rcpp::wrap(jr), 
                               Rcpp::Named("flux") = Rcpp::wrap(res_xr) );
};

Urban::Residents::Residents(Urban& urbs_, int from_): 
    urbs(urbs_), from(from_), col_dispo(map_col_dispo(urbs_, from_)) {};

Urban::SubRegion::SubRegion(Urban& urbs_, std::vector<int> gfrom, std::vector<int> gto):
    urbs(urbs_), group_from(gfrom), group_to(gto) {};

Urban::SubRegion::SubRegion(Urban& urbs_, Rcpp::IntegerVector gfrom, Rcpp::IntegerVector gto):
    urbs(urbs_), group_from(Rcpp::as< std::vector<int> >(gfrom)), group_to(Rcpp::as< std::vector<int> >(gto)) {};

std::vector<int> map_col_dispo(Urban& urbs_, int from_) {
    int debut = urbs_.p[from_];
    int fin = urbs_.p[from_ + 1];
    std::vector<int> col_dispo;
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
    int n_sites = res.col_dispo.size();
    std::vector<double> attract(n_sites);
    for (auto k = 0; k < n_sites; ++k) {
        auto index = res.col_dispo[k];
        attract[k] = res.urbs.emplois_libres[index] * (*fct)(res.urbs.xr[index]); // revoir def de la fonction
    }
    return attract;
}

std::vector<double> repartition_directe(Urban::Residents& res, 
                                  std::vector<double>& attract) {
  double total = std::accumulate(attract.begin(), attract.end(), 0.0);
  double absorption = -log(res.urbs.fuites[res.from]) / total;
  int n_sites = attract.size();
  
  // Calcul des actifs absorbés par sites.
  std::vector<double> jobtakers(attract);
  total = 0.0;
  double passants = res.urbs.actifs_libres[res.from];
  double anciens, facteur;
  
  for (auto index = 0; index < n_sites;) {
    auto pos = index + 1;
    while (pos != n_sites && res.urbs.xr[ res.col_dispo[index] ] == res.urbs.xr[ res.col_dispo[pos] ]) ++pos;
    total = std::accumulate(attract.begin() + index, attract.begin() + pos, 0.0);
    anciens = passants;
    passants *= exp(-absorption * total);
    facteur = (anciens - passants) / total;
    for (auto k = index; k < pos; ++k) {
      jobtakers[k] *= facteur;
    }
    index = pos;
  }
  return(jobtakers);
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
    int n_sites = cumul.size();

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

std::vector<float> Urban::SubRegion::regrouper(std::vector< std::vector<float> >& liens) {
    
    const int N = ng_from();
    const int K = ng_to();
    std::vector<float> flux(N * K);
  
#pragma omp declare reduction(vsumf : std::vector<float> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), \
    std::plus<float>())) initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp parallel for reduction(vsumf: flux) collapse(2)
    for (auto i = 0; i < urbs.n_from(); ++i) {
      for (auto j = 0; j < urbs.n_to(); ++j) {
        flux[ group_from[i] * K + group_to[j] ] += liens[i][j];
      }
    }
    return flux;
}

Rcpp::List Urban::SubRegion::procurer(std::vector< std::vector<float> >& liens) {

    const int N = ng_from();
    const int K = ng_to();
    const int taille = N * K;

    std::vector<float> res_xr(taille);
    res_xr = regrouper(liens);

    std::vector<int> res_i(taille);
    std::vector<int> res_j(taille);

    for (auto i = 0; i < N; ++i) {
        for (auto j = 0; j < K; ++j) {
            res_i[i * K + j] = i;
            res_j[i * K + j] = j;
      }
    }
    return Rcpp::List::create( Rcpp::Named("i") = Rcpp::wrap(res_i), 
                               Rcpp::Named("j") = Rcpp::wrap(res_j), 
                               Rcpp::Named("flux") = Rcpp::wrap(res_xr)); 
}

Rcpp::List Urban::SubRegion::procurer(std::vector< std::vector<float> >& liens,
                                     std::vector<float>& cible) {

    const int N = ng_from();
    const int K = ng_to();
    const int taille = N * K;

    std::vector<float> res_xr(taille);
    res_xr = regrouper(liens);

    std::vector<int> res_i(taille);
    std::vector<int> res_j(taille);

    for (auto i = 0; i < N; ++i) {
        for (auto j = 0; j < K; ++j) {
            res_i[i * K + j] = i;
            res_j[i * K + j] = j;
      }
    }

    float kl = entropie_relative<float>(res_xr, cible);

    return Rcpp::List::create( Rcpp::Named("i") = Rcpp::wrap(res_i), 
                               Rcpp::Named("j") = Rcpp::wrap(res_j), 
                               Rcpp::Named("flux") = Rcpp::wrap(res_xr),
                               Rcpp::Named("kl") = kl); 
}
