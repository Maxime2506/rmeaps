#ifndef __CLASSES__
#define __CLASSES__
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
    Urban();
    Urban(std::vector<int> jr, std::vector<int> p, std::vector<double> xr, 
        std::vector<double> actifs, std::vector<double> emplois, std::vector<double> fuites);

};

class Residents
{
private:
    const std::size_t from;
    std::vector<double> emplois_dispo;
    double actifs_libres;
    const double fuite;
    const Urban urb;
public:
    Residents();
    Residents(Urban urb, std::size_t from, std::vector<double> emplois_dispo, 
    double actifs_libres, double fuite);

};


template <class Iter>
class RowIter
{
private:
    Iter debut;
    Iter fin;
public:
    RowIter();
    RowIter(Iter d, Iter f);
    RowIter(const Urban region, const Iter line);

    Iter begin() { return debut;}
    Iter end() { return fin;}    
};

#endif // __CLASSES__