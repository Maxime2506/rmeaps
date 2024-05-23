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
};

std::vector<int> map_col_dispo(Urban& urbs_, int from_);
std::vector<double> attractivite(Urban::Residents& res, std::unique_ptr<fonction_attraction>);
std::vector<double> attract_cumul(Urban::Residents& res, std::vector<double>& attract);
std::vector<double> calc_repartition(Urban::Residents& res, std::vector<double>& attract, std::vector<double>& cumul);

#endif // __CLASSES__