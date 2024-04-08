#include <Rcpp.h>
//' @importClassFrom Matrix dgRMatrix
 //' 
 // Ranked Row Sparse Matrix Class
 class RankedRSMatrix {
  public:
   Rcpp::NumericVector xr;
   Rcpp::IntegerVector jr, p, dim;
   
   // constructors
   // This one takes the data ALREADY RANKED
   RankedRSMatrix(Rcpp::NumericVector xr, Rcpp::IntegerVector jr, Rcpp::IntegerVector p, Rcpp::IntegerVector dim)
     : xr(xr), jr(jr), p(p), dim(dim) {}
   RankedRSMatrix() = default;

 // This one takes a dgRMatrix class object (R pkg Matrix) NOT RANKED.
   RankedRSMatrix(Rcpp::S4 m) {
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
   int nrow() { return dim[0]; }
   int ncol() { return dim[1]; }
   int nvalid() { return xr.size(); };
   Rcpp::NumericVector& validvalues() { return xr; };
   Rcpp::IntegerVector& innerIndexPtr() { return jr; };
   Rcpp::IntegerVector& outerIndexPtr() { return p; };
   
   // element lookup at specific index
   double at(int row, int col) const {
     for (int i = p[row]; i < p[row + 1L]; ++i) {
       if (jr[i] == col) return xr[i];
     }
     return NA_REAL;
   }
   
   
   // Attention : les jr commencent à l'indice 0 même s'il s'agit de rang.
   Rcpp::S4 unrank() {
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
 };

// Exposition en R
RCPP_EXPOSED_AS(RankedRSMatrix);


RCPP_MODULE(RankedRowSparseMatrix) {
  
  Rcpp::class_<RankedRSMatrix>( "RankedRSMatrix" )
  .constructor()
  .constructor<Rcpp::NumericVector, Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::IntegerVector>("Construction directe à partir de lignes rangées.")
  .constructor<Rcpp::S4>("Construction à partir d'une Matrix::dgRMatrix.")
  .field("xr", &RankedRSMatrix::xr)
  .field("jr", &RankedRSMatrix::jr)
  .field("p", &RankedRSMatrix::p)
  .field("dim", &RankedRSMatrix::dim)
  .method("nrow", &RankedRSMatrix::nrow)
  .method("ncol", &RankedRSMatrix::ncol)
  .method("nvalid", &RankedRSMatrix::nvalid)
  .method("at", &RankedRSMatrix::at)
  .method("unrank", &RankedRSMatrix::unrank)
  ;
}
