#include <Rcpp.h>
#include "RankedRSMatrix.h"

RCPP_EXPOSED_CLASS(RankedRSMatrix)
  
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
  .method("rankby", &RankedRSMatrix::rankby)
  .method("unrank", &RankedRSMatrix::unrank)
  ;
}