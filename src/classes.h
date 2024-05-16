#ifndef __CLASSES__
#define __CLASSES__

#include<vector>
#include "classes_attraction.h"
class Urban
{
public:
    const std::vector<unsigned int>& jr; // Remaniant la convention des Row Sparse Matrix, innerindex des j rangés selon l'ordre des x à chaque ligne.
    const std::vector<std::size_t>& p; // outerstarts (convention RSM).
    const std::vector<double>& xr; // values (convention RSM).
    const std::vector<double>& actifs; // Tous les actifs au départ.
    const std::vector<double>& emplois; // Tous les emplois au départ.
    const std::vector<double>& fuites;
    std::vector<double> actifs_libres; // Les actifs non occupés à un moment de la distribution.
    std::vector<double> emplois_libres; // Les emplois encore libres à un moment de la distribution.
    Urban();
    Urban(std::vector<unsigned int>& jr_, std::vector<std::size_t>& p_, std::vector<double>& xr_, 
            std::vector<double>& actifs_, std::vector<double>& emplois_, std::vector<double>& fuites_);

    class Residents
    {
        public:
        const Urban& urbs;
        const unsigned int from;
        std::vector<unsigned int> col_dispo; // les indices des colonnes avec emplois dispo rangés par proximité pour ces résidents.

        Residents() = delete;
        Residents(Urban& urbs_, unsigned int from_);
        
        friend std::vector<unsigned int> map_col_dispo(Urban& urbs_, unsigned int from_);
        friend std::vector<double> attractivite(Urban::Residents& res, std::unique_ptr<fonction_attraction>);
        friend std::vector<double> attract_cumul(Urban::Residents& res, 
                                                 std::vector<double>& attract);
        friend std::vector<double> calc_repartition(Urban::Residents& res, 
                                                    std::vector<double>& attract,
                                                    std::vector<double>& cumul);
        
    };
};

#endif // __CLASSES__