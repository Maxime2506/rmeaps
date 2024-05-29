#ifndef __CLASSES__
#define __CLASSES__

#include <Rcpp.h>
#include <vector>
#include <cstddef>

#include "classes_attraction.h"

class Urban
{
public:
      const std::vector<int> jr; // Remaniant la convention des Row Sparse Matrix, innerindex des j rangés selon l'ordre des x à chaque ligne.
      const std::vector<int> p; // outerstarts (convention RSM).
      const std::vector<double> xr; // values (convention RSM).
      const std::vector<double> actifs; // Tous les actifs au départ.
      const std::vector<double> emplois; // Tous les emplois au départ.
      const std::vector<double> fuites;
     
      Urban();
      Urban(Rcpp::IntegerVector jr_, Rcpp::IntegerVector p_, Rcpp::NumericVector xr_, 
            Rcpp::NumericVector actifs_, Rcpp::NumericVector emplois_, Rcpp::NumericVector fuites_);
      Urban(const std::vector<int> jr_, const std::vector<int> p_, const std::vector<double> xr_, 
            const std::vector<double> actifs_, const std::vector<double> emplois_, const std::vector<double> fuites_);

      Rcpp::List format_sortie(std::vector< std::vector<float> >& liens);

      const int n_from() const { return actifs.size(); }
      const int n_to() const { return emplois.size(); }
      class Residents
      {
      public:
            const Urban& urb;
            const int from;
            
            Residents() = delete;
            Residents(Urban& urb_, int from_);

            std::vector<int> map_col_dispo(std::vector<double>& emplois_libres);
            std::vector<double> attractivite(std::vector<double>& emplois_libres, std::shared_ptr<fonction_attraction>& fct);
            std::vector<double> attractivite(std::vector<double>& emplois_libres, std::vector<int>& col_dispo, std::shared_ptr<fonction_attraction>& fct);
            std::vector<double> repartition_nolimit(std::vector<int>& col_dispo, std::vector<double>& attract, double actifs_restants);
            std::vector<double> repartition_limited(std::vector<double>& emplois_libres, std::vector<int> col_dispo, std::vector<double> attract);
      };

      class SubRegion
      {
      public:
            const Urban& urb;
            const std::vector<int> group_from, group_to;

            SubRegion() = delete;
            SubRegion(Urban& urb_, std::vector<int> gfrom, std::vector<int> gto);
            SubRegion(Urban& urb_, Rcpp::IntegerVector gfrom, Rcpp::IntegerVector gto);

            const int ng_from() const { return 1L + *std::max_element(group_from.begin(), group_from.end()); };
            const int ng_to() const { return 1L + *std::max_element(group_to.begin(), group_to.end()); };

            std::vector<float> regrouper(std::vector< std::vector<float> >& liens);
            Rcpp::List format_sortie(std::vector< std::vector<float> >& liens);
            Rcpp::List format_sortie(std::vector< std::vector<float> >& liens, std::vector<float>& cible);
      };
};

#endif // __CLASSES__