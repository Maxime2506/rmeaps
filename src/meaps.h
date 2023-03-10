#ifndef __MEAPS__
#define __MEAPS__

Rcpp::NumericMatrix meaps_oneshuf(
    const Rcpp::IntegerMatrix rkdist, 
    Rcpp::NumericVector emplois,
    const Rcpp::NumericVector actifs,
    const Rcpp::NumericMatrix modds,
    const Rcpp::NumericVector f, 
    const Rcpp::IntegerVector shuf,
    bool normalisation,
    double fuite_min,
    double seuil_newton);

#endif // __MEAPS__