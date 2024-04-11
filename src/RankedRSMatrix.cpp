#include <Rcpp.h>
#include "RankedRSMatrix.h"

//' @importClassFrom Matrix dgRMatrix
//' 
// Ranked Row Sparse Matrix Class
// class RankedRSMatrix 
 
// constructors
RankedRSMatrix::RankedRSMatrix() = default;
// This one takes the data ALREADY RANKED
RankedRSMatrix::RankedRSMatrix(Rcpp::NumericVector xr, Rcpp::IntegerVector jr, Rcpp::IntegerVector p, Rcpp::IntegerVector dim)
     : xr(xr), jr(jr), p(p), dim(dim) {}

// This one takes a dgRMatrix class object (R pkg Matrix) NOT RANKED.
RankedRSMatrix::RankedRSMatrix(Rcpp::S4 m) {
     Rcpp::NumericVector xr_ = m.slot("x");
     Rcpp::IntegerVector jr_ = m.slot("j");
     Rcpp::IntegerVector p_ = m.slot("p");
     Rcpp::IntegerVector dim_ = m.slot("Dim");
     
     dim = Rcpp::clone(dim_);
     p = Rcpp::clone(p_);
     jr = Rcpp::clone(jr_);
     xr = Rcpp::clone( xr_ );
     
     for (auto from = 0; from < dim[0]; ++from) {
       std::multimap<double, int> arrangement;
       int debut = p[from], fin = p[from + 1L];
       for (auto k = debut; k < fin; ++k) {
         arrangement.insert(std::make_pair(xr(k), jr(k)));
       }
       
       int index = 0L;
       for (auto it = arrangement.begin(); it != arrangement.end(); ++it) {
         xr[debut + index] = (it->first);
         jr[debut + index] = (it->second);
         ++index;
       }
     }
   }
   
// methods
int RankedRSMatrix::nrow() { return dim[0]; }
int RankedRSMatrix::ncol() { return dim[1]; }
int RankedRSMatrix::nvalid() { return xr.size(); };
Rcpp::NumericVector& RankedRSMatrix::validvalues() { return xr; };
Rcpp::IntegerVector& RankedRSMatrix::innerIndexPtr() { return jr; };
Rcpp::IntegerVector& RankedRSMatrix::outerIndexPtr() { return p; };
   
// element lookup at specific index
double RankedRSMatrix::at(int row, int col) {
  for (int i = p[row]; i < p[row + 1L]; ++i) {
    if (jr[i] == col) return xr[i];
  }
  return NA_REAL;
};
   
// method pour changer l'ordre jr selon les rangs d'une autre sparse matrix
void RankedRSMatrix::rankby(RankedRSMatrix d) {
  
  std::vector<int> new_jr;
  std::vector<double> new_xr;
  
  for (auto from = 0; from < dim[0]; ++from) {
    
    std::map<int, int> rangementby;
    int debut = d.p[from], fin = d.p[from + 1L];
    for (auto k = p[from]; k < p[from + 1L]; ++k) {
      auto pos = std::find(d.jr.begin() + debut, d.jr.begin() + fin, jr[k]);
      if (pos != d.jr.begin() + fin) {
        auto ind = std::distance(d.jr.begin() + debut, pos);
        rangementby.insert(std::make_pair(ind, k));
      } else {
        Rcpp::warning("Des odds valides renvoient vers des distances invalides.");
      }
    }
    
    for (auto it = rangementby.begin(); it != rangementby.end(); ++it) {
      std::size_t k = (it->second);
      new_jr.push_back(jr[k]);
      new_xr.push_back(xr[k]);
    }
  }
  jr = Rcpp::wrap(new_jr);
  xr = Rcpp::wrap(new_xr);
}
   
// Attention : les jr commencent à l'indice 0 même s'il s'agit de rang.
Rcpp::S4 RankedRSMatrix::unrank() {
  Rcpp::NumericVector x(xr.size());
  std::vector<int> j = Rcpp::as< std::vector<int> >(Rcpp::clone(jr));
  
  for (auto from = 0; from < dim[0]; ++from) {
    int debut = p[from], fin = p[from + 1L];
    std::sort(j.begin() + debut, j.begin() + fin);
    
    for (auto k = debut; k < fin; ++k) {
      auto pos = std::lower_bound(j.begin() + debut, j.begin() + fin, jr[k]);
      auto ind = std::distance(j.begin(), pos);
      x[ind] = xr[k];
    }
  }
  
  Rcpp::S4 res (std::string("dgRMatrix"));
  res.slot("p") = Rcpp::clone(p);
  res.slot("j") = Rcpp::wrap(j);
  res.slot("Dim") = Rcpp::clone(dim);
  res.slot("x") = Rcpp::wrap(x);
  
  return res;
};



