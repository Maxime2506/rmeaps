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
        const unsigned int from;
        const std::vector<double>& emplois_dispo;
        const double actifs_libres;
    public:
        Residents() = default;
        Residents(Urban& urbs_, std::size_t from_, std::vector<double>& emplois_dispo_, double actifs_libres_);
        Residents(Urban& urbs_, std::size_t from_);
    };
    friend class InnerIterator;
    friend class StrictInnerIterator;
};

Urban::Urban(std::vector<int>& jr_, std::vector<int>& p_, std::vector<double>& xr_, 
        std::vector<double>& actifs_, std::vector<double>& emplois_, std::vector<double>& fuites_):
        jr(jr_), p(p_), xr(xr_), actifs(actifs_), emplois(emplois_), fuites(fuites_) {};

Urban::Residents::Residents(Urban& urbs_, std::size_t from_, std::vector<double>& emplois_dispo_, double actifs_libres_):
        urbs(urbs_), from(from_), emplois_dispo(emplois_dispo_), actifs_libres(actifs_libres_) {};

Urban::Residents::Residents(Urban& urbs_, std::size_t from_):
        urbs(urbs_), from(from_), emplois_dispo(urbs_.emplois), actifs_libres(urbs_.actifs[from_]) {};


class InnerIterator {
private:
    const Urban& ptr;
    const unsigned int from, max_index;
    unsigned int index;
public:
    InnerIterator(Urban& ptr, unsigned int row) : ptr(ptr), from(row), index(ptr.p[row]), max_index(ptr.p[row + 1]) {}
    operator bool() const { return (index < max_index); }
    InnerIterator& operator++() {
        ++index;
        return *this;
    }
    const double& value() const { return ptr.xr[index]; }
    unsigned int col() const { return ptr.jr[index]; }
    unsigned int row() const { return from; }
};


class StrictInnerIterator {
private:
    const Urban& ptr;
    const unsigned int from, max_index;
    unsigned int index;
public:
    StrictInnerIterator(Urban& ptr, unsigned int row) : ptr(ptr), from(row), index(ptr.p[row]), max_index(ptr.p[row + 1]) {}
    operator bool() const { return (index < max_index); }
    StrictInnerIterator& operator++() {
        while(index < max_index && ptr.xr[index] == ptr.xr[index + 1]) ++index;
        return *this;
    }
    const double cumsum(const std::vector<double>& values, const StrictInnerIterator& pos) const { 
        std::size_t fin = pos.index;
        while (fin < pos.max_index && pos.ptr.xr[fin] == pos.ptr.xr[fin + 1]) ++fin;
        double tot = 0.0;
        for (auto equal_rank = pos.index; equal_rank < fin; ++equal_rank) tot += values[equal_rank];
        return tot; 
        }
};

