#ifndef __RankedRSMatrix__
#define __RankedRSMatrix__

class RankedRSMatrix {
public:
  Rcpp::NumericVector xr;
  Rcpp::IntegerVector jr, p, dim;
  
  RankedRSMatrix();
  RankedRSMatrix(Rcpp::NumericVector xr, Rcpp::IntegerVector jr, Rcpp::IntegerVector p, Rcpp::IntegerVector dim);
  RankedRSMatrix(Rcpp::S4 m);
  
  int nrow();
  int ncol();
  int nvalid();
  Rcpp::NumericVector& validvalues();
  Rcpp::IntegerVector& innerIndexPtr();
  Rcpp::IntegerVector& outerIndexPtr();
  
  double at(int row, int col);

  void rankby(RankedRSMatrix d);
  Rcpp::S4 unrank();
};

#endif // __RankedRSMatrix__
