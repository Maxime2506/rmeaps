// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// deborder
List deborder(NumericVector conteneurs, const double& quantite);
RcppExport SEXP _rmeaps_deborder(SEXP conteneursSEXP, SEXP quantiteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type conteneurs(conteneursSEXP);
    Rcpp::traits::input_parameter< const double& >::type quantite(quantiteSEXP);
    rcpp_result_gen = Rcpp::wrap(deborder(conteneurs, quantite));
    return rcpp_result_gen;
END_RCPP
}
// distribuer
List distribuer(NumericVector conteneurs, NumericVector proportion, const double& quantite);
RcppExport SEXP _rmeaps_distribuer(SEXP conteneursSEXP, SEXP proportionSEXP, SEXP quantiteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type conteneurs(conteneursSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type proportion(proportionSEXP);
    Rcpp::traits::input_parameter< const double& >::type quantite(quantiteSEXP);
    rcpp_result_gen = Rcpp::wrap(distribuer(conteneurs, proportion, quantite));
    return rcpp_result_gen;
END_RCPP
}
// meaps_single
NumericMatrix meaps_single(const IntegerMatrix rkdist, const NumericVector emplois, const NumericVector actifs, const NumericMatrix modds, const NumericVector f, IntegerVector shuf);
RcppExport SEXP _rmeaps_meaps_single(SEXP rkdistSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP moddsSEXP, SEXP fSEXP, SEXP shufSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type rkdist(rkdistSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type modds(moddsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type shuf(shufSEXP);
    rcpp_result_gen = Rcpp::wrap(meaps_single(rkdist, emplois, actifs, modds, f, shuf));
    return rcpp_result_gen;
END_RCPP
}
// meaps_bootstrap2
NumericMatrix meaps_bootstrap2(IntegerMatrix rkdist, NumericVector emplois, NumericVector actifs, NumericMatrix modds, NumericVector f, IntegerMatrix shuf);
RcppExport SEXP _rmeaps_meaps_bootstrap2(SEXP rkdistSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP moddsSEXP, SEXP fSEXP, SEXP shufSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type rkdist(rkdistSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type modds(moddsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type shuf(shufSEXP);
    rcpp_result_gen = Rcpp::wrap(meaps_bootstrap2(rkdist, emplois, actifs, modds, f, shuf));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rmeaps_deborder", (DL_FUNC) &_rmeaps_deborder, 2},
    {"_rmeaps_distribuer", (DL_FUNC) &_rmeaps_distribuer, 3},
    {"_rmeaps_meaps_single", (DL_FUNC) &_rmeaps_meaps_single, 6},
    {"_rmeaps_meaps_bootstrap2", (DL_FUNC) &_rmeaps_meaps_bootstrap2, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_rmeaps(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
