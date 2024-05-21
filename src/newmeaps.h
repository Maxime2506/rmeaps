#ifndef __MEAPS_CLASS__
#define __MEAPS_CLASS__

Rcpp::List meapsclass(const Rcpp::IntegerVector jr_dist, 
                const Rcpp::IntegerVector p_dist, 
                const Rcpp::NumericVector xr_dist, 
                const Rcpp::NumericVector emplois,
                const Rcpp::NumericVector actifs, 
                const Rcpp::NumericVector fuites, 
                const Rcpp::NumericVector parametres,
                const std::string attraction = "constant",
                const Rcpp::Nullable<Rcpp::IntegerVector> group_from = R_NilValue,
                const Rcpp::Nullable<Rcpp::IntegerVector> group_to = R_NilValue,
                const Rcpp::Nullable<Rcpp::NumericVector> cible = R_NilValue,
                const int nthreads = 0L, 
                const bool verbose = true);

#endif // __MEAPS_CLASS__