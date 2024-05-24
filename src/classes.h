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
      std::vector<double> actifs_libres; // Les actifs non occupés à un moment de la distribution.
      std::vector<double> emplois_libres; // Les emplois encore libres à un moment de la distribution.
      Urban();
      Urban(Rcpp::IntegerVector jr_, Rcpp::IntegerVector p_, Rcpp::NumericVector xr_, 
            Rcpp::NumericVector actifs_, Rcpp::NumericVector emplois_, Rcpp::NumericVector fuites_);
      Urban(const std::vector<int> jr_, const std::vector<int> p_, const std::vector<double> xr_, 
            const std::vector<double> actifs_, const std::vector<double> emplois_, const std::vector<double> fuites_);

      Rcpp::List procurer(std::vector< std::vector<float> >& liens);

      const int n_from() const { return actifs.size(); }
      const int n_to() const { return emplois.size(); }
      class Residents
      {
      public:
            const Urban& urbs;
            const int from;
            std::vector<int> col_dispo; // les indices des colonnes avec emplois dispo rangés par proximité pour ces résidents.

            Residents() = delete;
            Residents(Urban& urbs_, int from_);
        
            friend std::vector<int> map_col_dispo(Urban& urbs_, int from_);
            friend std::vector<double> attractivite(Urban::Residents& res, std::shared_ptr<fonction_attraction>& fct);
            friend std::vector<double> repartition_directe(Urban::Residents& res, 
                                                 std::vector<double>& attract);
            friend std::vector<double> attract_cumul(Urban::Residents& res, 
                                                 std::vector<double>& attract);
            friend std::vector<double> calc_repartition(Urban::Residents& res, 
                                                    std::vector<double>& attract,
                                                    std::vector<double>& cumul);
        
      };

      class SubRegion
      {
      public:
            const Urban& urbs;
            const std::vector<int> group_from, group_to;

            SubRegion() = delete;
            SubRegion(Urban& urbs_, std::vector<int> gfrom, std::vector<int> gto);
            SubRegion(Urban& urbs_, Rcpp::IntegerVector gfrom, Rcpp::IntegerVector gto);

            const int ng_from() const { return 1L + *std::max_element(group_from.begin(), group_from.end()); };
            const int ng_to() const { return 1L + *std::max_element(group_to.begin(), group_to.end()); };

            std::vector<float> regrouper(std::vector< std::vector<float> >& liens);
            Rcpp::List procurer(std::vector< std::vector<float> >& liens);
            Rcpp::List procurer(std::vector< std::vector<float> >& liens, std::vector<float>& cible);
      };
};

std::vector<int> map_col_dispo(Urban& urbs_, int from_);
std::vector<double> attractivite(Urban::Residents& res, std::unique_ptr<fonction_attraction>& fct);
std::vector<double> repartition_directe(Urban::Residents& res, std::vector<double>& attract);
std::vector<double> attract_cumul(Urban::Residents& res, std::vector<double>& attract);
std::vector<double> calc_repartition(Urban::Residents& res, std::vector<double>& attract, std::vector<double>& cumul);

#endif // __CLASSES__