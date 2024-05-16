#include <type_traits>
#include <vector>
#include <iterator>

#include <iostream>
using namespace std;

class Urban
{
private:
    const std::vector<int>& jr; // Remaniant la convention des Row Sparse Matrix, innerindex des j rangés selon l'ordre des x à chaque ligne.
    const std::vector<int>& p; // outerstarts (convention RSM).
    const std::vector<double>& xr; // values (convention RSM)
    const std::vector<double>& actifs; 
    const std::vector<double>& emplois;
    const std::vector<double>& fuites;
public:
    Urban() = default;
    Urban(std::vector<int>& jr_, std::vector<int>& p_, std::vector<double>& xr_, 
            std::vector<double>& actifs_, std::vector<double>& emplois_, std::vector<double>& fuites_);

    class Residents
    {
    private:
        const Urban& urbs;
        const std::size_t from;
        const std::vector<double>& emplois_dispo;
        const double actifs_libres;
    public:
        Residents() = default;
        Residents(Urban& urbs_, std::size_t from_, std::vector<double>& emplois_dispo_, double actifs_libres_);
        Residents(Urban& urbs_, std::size_t from_);
    };
    //friend RowIter<class Line::iterator> parcourir(const Urban& region, const std::size_t line);
};

Urban::Urban(std::vector<int>& jr_, std::vector<int>& p_, std::vector<double>& xr_, 
        std::vector<double>& actifs_, std::vector<double>& emplois_, std::vector<double>& fuites_):
        jr(jr_), p(p_), xr(xr_), actifs(actifs_), emplois(emplois_), fuites(fuites_) {};

Urban::Residents::Residents(Urban& urbs_, std::size_t from_, std::vector<double>& emplois_dispo_, double actifs_libres_):
        urbs(urbs_), from(from_), emplois_dispo(emplois_dispo_), actifs_libres(actifs_libres_) {};
Urban::Residents::Residents(Urban& urbs_, std::size_t from_):
        urbs(urbs_), from(from_), emplois_dispo(urbs_.emplois), actifs_libres(urbs_.actifs[from_]) {};

template <typename Iter> 
class RowIter
{
    Iter debut;
    Iter fin;
public:
    RowIter() = default;
    RowIter(Iter d, Iter f);
    RowIter(const Urban region, const Iter line);

    Iter begin() { return debut;}
    Iter end() { return fin;}    
};

template <typename Iter>
RowIter<Iter>::RowIter(Iter d, Iter f): debut(d), fin(f) {};

template <typename Iter>
RowIter<Iter>::RowIter(const Urban region, const Iter line)
{
    debut = region.p[line];
    fin = region.p[line + 1L];
}

template <typename Line>
RowIter<typename Line::iterator> parcourir(const Urban& region, typename line)
{
    return RowIter<typename line::iterator> (c.begin() + debut, c.begin() + fin);
}


