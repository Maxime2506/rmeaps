#ifndef __MEAPS__
#define __MEAPS__

Rcpp::NumericMatrix meaps_single(
    const Rcpp::IntegerMatrix rkdist, 
    const Rcpp::NumericVector emplois,
    const Rcpp::NumericVector actifs,
    const Rcpp::NumericMatrix modds,
    const Rcpp::NumericVector f, 
    Rcpp::IntegerVector shuf);

#endif // __MEAPS__