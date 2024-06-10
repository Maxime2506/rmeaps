// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// meaps_all_in_cpp
List meaps_all_in_cpp(const IntegerVector jr_dist, const IntegerVector p_dist, const NumericVector xr_dist, const NumericVector emplois, const NumericVector actifs, const NumericVector fuites, const NumericVector parametres, const std::string attraction, const Nullable<IntegerVector> group_from, const Nullable<IntegerVector> group_to, const Nullable<NumericVector> cible, const int nthreads, const bool verbose);
RcppExport SEXP _rmeaps_meaps_all_in_cpp(SEXP jr_distSEXP, SEXP p_distSEXP, SEXP xr_distSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP fuitesSEXP, SEXP parametresSEXP, SEXP attractionSEXP, SEXP group_fromSEXP, SEXP group_toSEXP, SEXP cibleSEXP, SEXP nthreadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type jr_dist(jr_distSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type p_dist(p_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type xr_dist(xr_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type fuites(fuitesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type parametres(parametresSEXP);
    Rcpp::traits::input_parameter< const std::string >::type attraction(attractionSEXP);
    Rcpp::traits::input_parameter< const Nullable<IntegerVector> >::type group_from(group_fromSEXP);
    Rcpp::traits::input_parameter< const Nullable<IntegerVector> >::type group_to(group_toSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector> >::type cible(cibleSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(meaps_all_in_cpp(jr_dist, p_dist, xr_dist, emplois, actifs, fuites, parametres, attraction, group_from, group_to, cible, nthreads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// meapsmode_cpp
List meapsmode_cpp(const IntegerVector jr_dist, const IntegerVector p_dist, const NumericVector xr_dist, const NumericVector emplois, const NumericVector actifs, const NumericVector fuites, const IntegerVector j_mode, const NumericVector parametres, const std::string attraction, const Nullable<IntegerVector> group_from, const Nullable<IntegerVector> group_to, const Nullable<NumericVector> cible, const int nthreads, const bool verbose);
RcppExport SEXP _rmeaps_meapsmode_cpp(SEXP jr_distSEXP, SEXP p_distSEXP, SEXP xr_distSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP fuitesSEXP, SEXP j_modeSEXP, SEXP parametresSEXP, SEXP attractionSEXP, SEXP group_fromSEXP, SEXP group_toSEXP, SEXP cibleSEXP, SEXP nthreadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type jr_dist(jr_distSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type p_dist(p_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type xr_dist(xr_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type fuites(fuitesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type j_mode(j_modeSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type parametres(parametresSEXP);
    Rcpp::traits::input_parameter< const std::string >::type attraction(attractionSEXP);
    Rcpp::traits::input_parameter< const Nullable<IntegerVector> >::type group_from(group_fromSEXP);
    Rcpp::traits::input_parameter< const Nullable<IntegerVector> >::type group_to(group_toSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector> >::type cible(cibleSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(meapsmode_cpp(jr_dist, p_dist, xr_dist, emplois, actifs, fuites, j_mode, parametres, attraction, group_from, group_to, cible, nthreads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// multishuf_oc_group_cpp
List multishuf_oc_group_cpp(const IntegerVector jr_dist, const IntegerVector p_dist, const NumericVector xr_dist, const NumericVector emplois, const NumericVector actifs, const NumericVector fuites, const IntegerMatrix shuf, const IntegerVector group_from, const IntegerVector group_to, const NumericVector parametres, const NumericVector xr_odds, const Nullable<NumericVector> cible, const std::string attraction, int nthreads, bool verbose);
RcppExport SEXP _rmeaps_multishuf_oc_group_cpp(SEXP jr_distSEXP, SEXP p_distSEXP, SEXP xr_distSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP fuitesSEXP, SEXP shufSEXP, SEXP group_fromSEXP, SEXP group_toSEXP, SEXP parametresSEXP, SEXP xr_oddsSEXP, SEXP cibleSEXP, SEXP attractionSEXP, SEXP nthreadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type jr_dist(jr_distSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type p_dist(p_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type xr_dist(xr_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type fuites(fuitesSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type shuf(shufSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type group_from(group_fromSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type group_to(group_toSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type parametres(parametresSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type xr_odds(xr_oddsSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector> >::type cible(cibleSEXP);
    Rcpp::traits::input_parameter< const std::string >::type attraction(attractionSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(multishuf_oc_group_cpp(jr_dist, p_dist, xr_dist, emplois, actifs, fuites, shuf, group_from, group_to, parametres, xr_odds, cible, attraction, nthreads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// max_threads
IntegerVector max_threads();
RcppExport SEXP _rmeaps_max_threads() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(max_threads());
    return rcpp_result_gen;
END_RCPP
}
// multishuf_oc_cpp
List multishuf_oc_cpp(const IntegerVector jr_dist, const IntegerVector p_dist, const NumericVector xr_dist, const NumericVector emplois, const NumericVector actifs, const NumericVector fuites, const IntegerMatrix shuf, const NumericVector parametres, const NumericVector xr_odds, const std::string attraction, int nthreads, bool verbose);
RcppExport SEXP _rmeaps_multishuf_oc_cpp(SEXP jr_distSEXP, SEXP p_distSEXP, SEXP xr_distSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP fuitesSEXP, SEXP shufSEXP, SEXP parametresSEXP, SEXP xr_oddsSEXP, SEXP attractionSEXP, SEXP nthreadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type jr_dist(jr_distSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type p_dist(p_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type xr_dist(xr_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type fuites(fuitesSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type shuf(shufSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type parametres(parametresSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type xr_odds(xr_oddsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type attraction(attractionSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(multishuf_oc_cpp(jr_dist, p_dist, xr_dist, emplois, actifs, fuites, shuf, parametres, xr_odds, attraction, nthreads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// multishuf_origin_cpp
List multishuf_origin_cpp(const IntegerVector jr_dist, const IntegerVector p_dist, const NumericVector xr_dist, const NumericVector emplois, const NumericVector actifs, const NumericVector fuites, const IntegerMatrix shuf, const std::string attraction, const NumericVector parametres, const NumericVector xr_odds, const std::string mode, const Nullable<NumericVector> oddssubjectifs, const int nthreads, const bool verbose);
RcppExport SEXP _rmeaps_multishuf_origin_cpp(SEXP jr_distSEXP, SEXP p_distSEXP, SEXP xr_distSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP fuitesSEXP, SEXP shufSEXP, SEXP attractionSEXP, SEXP parametresSEXP, SEXP xr_oddsSEXP, SEXP modeSEXP, SEXP oddssubjectifsSEXP, SEXP nthreadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type jr_dist(jr_distSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type p_dist(p_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type xr_dist(xr_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type fuites(fuitesSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type shuf(shufSEXP);
    Rcpp::traits::input_parameter< const std::string >::type attraction(attractionSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type parametres(parametresSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type xr_odds(xr_oddsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector> >::type oddssubjectifs(oddssubjectifsSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(multishuf_origin_cpp(jr_dist, p_dist, xr_dist, emplois, actifs, fuites, shuf, attraction, parametres, xr_odds, mode, oddssubjectifs, nthreads, verbose));
    return rcpp_result_gen;
END_RCPP
}
// multishuf_task_cpp
List multishuf_task_cpp(const IntegerVector jr_dist, const IntegerVector p_dist, const NumericVector xr_dist, const NumericVector emplois, const NumericVector actifs, const NumericVector fuites, const NumericVector parametres, const IntegerMatrix shuf, const std::string attraction, const Nullable<IntegerVector> group_from, const Nullable<IntegerVector> group_to, const Nullable<NumericVector> cible, const int nthreads, const bool verbose);
RcppExport SEXP _rmeaps_multishuf_task_cpp(SEXP jr_distSEXP, SEXP p_distSEXP, SEXP xr_distSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP fuitesSEXP, SEXP parametresSEXP, SEXP shufSEXP, SEXP attractionSEXP, SEXP group_fromSEXP, SEXP group_toSEXP, SEXP cibleSEXP, SEXP nthreadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type jr_dist(jr_distSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type p_dist(p_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type xr_dist(xr_distSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type fuites(fuitesSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type parametres(parametresSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type shuf(shufSEXP);
    Rcpp::traits::input_parameter< const std::string >::type attraction(attractionSEXP);
    Rcpp::traits::input_parameter< const Nullable<IntegerVector> >::type group_from(group_fromSEXP);
    Rcpp::traits::input_parameter< const Nullable<IntegerVector> >::type group_to(group_toSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector> >::type cible(cibleSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(multishuf_task_cpp(jr_dist, p_dist, xr_dist, emplois, actifs, fuites, parametres, shuf, attraction, group_from, group_to, cible, nthreads, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rmeaps_meaps_all_in_cpp", (DL_FUNC) &_rmeaps_meaps_all_in_cpp, 13},
    {"_rmeaps_meapsmode_cpp", (DL_FUNC) &_rmeaps_meapsmode_cpp, 14},
    {"_rmeaps_multishuf_oc_group_cpp", (DL_FUNC) &_rmeaps_multishuf_oc_group_cpp, 15},
    {"_rmeaps_max_threads", (DL_FUNC) &_rmeaps_max_threads, 0},
    {"_rmeaps_multishuf_oc_cpp", (DL_FUNC) &_rmeaps_multishuf_oc_cpp, 12},
    {"_rmeaps_multishuf_origin_cpp", (DL_FUNC) &_rmeaps_multishuf_origin_cpp, 14},
    {"_rmeaps_multishuf_task_cpp", (DL_FUNC) &_rmeaps_multishuf_task_cpp, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_rmeaps(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
