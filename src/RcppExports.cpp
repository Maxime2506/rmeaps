// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// chances_absorption
NumericMatrix chances_absorption(const IntegerMatrix rkdist, const NumericVector emplois, const NumericMatrix modds, const NumericVector f);
RcppExport SEXP _rmeaps_chances_absorption(SEXP rkdistSEXP, SEXP emploisSEXP, SEXP moddsSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type rkdist(rkdistSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type modds(moddsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(chances_absorption(rkdist, emplois, modds, f));
    return rcpp_result_gen;
END_RCPP
}
// communaliser
NumericMatrix communaliser(NumericMatrix flux, IntegerVector group_orig, IntegerVector group_dest);
RcppExport SEXP _rmeaps_communaliser(SEXP fluxSEXP, SEXP group_origSEXP, SEXP group_destSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type flux(fluxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group_orig(group_origSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group_dest(group_destSEXP);
    rcpp_result_gen = Rcpp::wrap(communaliser(flux, group_orig, group_dest));
    return rcpp_result_gen;
END_RCPP
}
// meaps_oneshuf
NumericMatrix meaps_oneshuf(IntegerMatrix rkdist, NumericVector emplois, NumericVector actifs, NumericMatrix modds, NumericVector f, IntegerVector shuf, std::string mode, Nullable<NumericVector> oddssubjectifs, bool normalisation, double fuite_min, double seuil_newton);
RcppExport SEXP _rmeaps_meaps_oneshuf(SEXP rkdistSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP moddsSEXP, SEXP fSEXP, SEXP shufSEXP, SEXP modeSEXP, SEXP oddssubjectifsSEXP, SEXP normalisationSEXP, SEXP fuite_minSEXP, SEXP seuil_newtonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type rkdist(rkdistSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type modds(moddsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type shuf(shufSEXP);
    Rcpp::traits::input_parameter< std::string >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type oddssubjectifs(oddssubjectifsSEXP);
    Rcpp::traits::input_parameter< bool >::type normalisation(normalisationSEXP);
    Rcpp::traits::input_parameter< double >::type fuite_min(fuite_minSEXP);
    Rcpp::traits::input_parameter< double >::type seuil_newton(seuil_newtonSEXP);
    rcpp_result_gen = Rcpp::wrap(meaps_oneshuf(rkdist, emplois, actifs, modds, f, shuf, mode, oddssubjectifs, normalisation, fuite_min, seuil_newton));
    return rcpp_result_gen;
END_RCPP
}
// meaps_continu_cpp
NumericVector meaps_continu_cpp(IntegerVector j_dist, IntegerVector p_dist, NumericVector x_dist, NumericVector emplois, NumericVector actifs, NumericVector f, IntegerMatrix shuf, NumericVector param, NumericVector j_odds, NumericVector p_odds, NumericVector x_odds, std::string attraction, int nthreads, bool progress, bool normalisation, double fuite_min);
RcppExport SEXP _rmeaps_meaps_continu_cpp(SEXP j_distSEXP, SEXP p_distSEXP, SEXP x_distSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP fSEXP, SEXP shufSEXP, SEXP paramSEXP, SEXP j_oddsSEXP, SEXP p_oddsSEXP, SEXP x_oddsSEXP, SEXP attractionSEXP, SEXP nthreadsSEXP, SEXP progressSEXP, SEXP normalisationSEXP, SEXP fuite_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type j_dist(j_distSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p_dist(p_distSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_dist(x_distSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type shuf(shufSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type j_odds(j_oddsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p_odds(p_oddsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_odds(x_oddsSEXP);
    Rcpp::traits::input_parameter< std::string >::type attraction(attractionSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< bool >::type normalisation(normalisationSEXP);
    Rcpp::traits::input_parameter< double >::type fuite_min(fuite_minSEXP);
    rcpp_result_gen = Rcpp::wrap(meaps_continu_cpp(j_dist, p_dist, x_dist, emplois, actifs, f, shuf, param, j_odds, p_odds, x_odds, attraction, nthreads, progress, normalisation, fuite_min));
    return rcpp_result_gen;
END_RCPP
}
// meaps_multishuf
NumericMatrix meaps_multishuf(IntegerMatrix rkdist, NumericVector emplois, NumericVector actifs, NumericMatrix modds, NumericVector f, IntegerMatrix shuf, std::string mode, Nullable<NumericVector> oddssubjectifs, int nthreads, bool progress, bool normalisation, double fuite_min, double seuil_newton);
RcppExport SEXP _rmeaps_meaps_multishuf(SEXP rkdistSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP moddsSEXP, SEXP fSEXP, SEXP shufSEXP, SEXP modeSEXP, SEXP oddssubjectifsSEXP, SEXP nthreadsSEXP, SEXP progressSEXP, SEXP normalisationSEXP, SEXP fuite_minSEXP, SEXP seuil_newtonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type rkdist(rkdistSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type modds(moddsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type shuf(shufSEXP);
    Rcpp::traits::input_parameter< std::string >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type oddssubjectifs(oddssubjectifsSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< bool >::type normalisation(normalisationSEXP);
    Rcpp::traits::input_parameter< double >::type fuite_min(fuite_minSEXP);
    Rcpp::traits::input_parameter< double >::type seuil_newton(seuil_newtonSEXP);
    rcpp_result_gen = Rcpp::wrap(meaps_multishuf(rkdist, emplois, actifs, modds, f, shuf, mode, oddssubjectifs, nthreads, progress, normalisation, fuite_min, seuil_newton));
    return rcpp_result_gen;
END_RCPP
}
// meaps_optim_cpp
NumericMatrix meaps_optim_cpp(IntegerVector jr_dist, IntegerVector p_dist, NumericVector xr_dist, NumericVector emplois, NumericVector actifs, NumericVector f, IntegerMatrix shuf, IntegerVector row_group, IntegerVector col_group, NumericVector param, IntegerVector jr_odds, IntegerVector p_odds, NumericVector xr_odds, std::string attraction, int nthreads, bool progress, bool normalisation, double fuite_min);
RcppExport SEXP _rmeaps_meaps_optim_cpp(SEXP jr_distSEXP, SEXP p_distSEXP, SEXP xr_distSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP fSEXP, SEXP shufSEXP, SEXP row_groupSEXP, SEXP col_groupSEXP, SEXP paramSEXP, SEXP jr_oddsSEXP, SEXP p_oddsSEXP, SEXP xr_oddsSEXP, SEXP attractionSEXP, SEXP nthreadsSEXP, SEXP progressSEXP, SEXP normalisationSEXP, SEXP fuite_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type jr_dist(jr_distSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p_dist(p_distSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xr_dist(xr_distSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type shuf(shufSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type row_group(row_groupSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type col_group(col_groupSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type jr_odds(jr_oddsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p_odds(p_oddsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xr_odds(xr_oddsSEXP);
    Rcpp::traits::input_parameter< std::string >::type attraction(attractionSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< bool >::type normalisation(normalisationSEXP);
    Rcpp::traits::input_parameter< double >::type fuite_min(fuite_minSEXP);
    rcpp_result_gen = Rcpp::wrap(meaps_optim_cpp(jr_dist, p_dist, xr_dist, emplois, actifs, f, shuf, row_group, col_group, param, jr_odds, p_odds, xr_odds, attraction, nthreads, progress, normalisation, fuite_min));
    return rcpp_result_gen;
END_RCPP
}
// one_distrib_continu
std::vector<double> one_distrib_continu(const double& entrants, const double& fuite, const std::vector<double>& attractivite, const std::vector<double>& distances, const std::vector<double>& placeslibres);
RcppExport SEXP _rmeaps_one_distrib_continu(SEXP entrantsSEXP, SEXP fuiteSEXP, SEXP attractiviteSEXP, SEXP distancesSEXP, SEXP placeslibresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type entrants(entrantsSEXP);
    Rcpp::traits::input_parameter< const double& >::type fuite(fuiteSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type attractivite(attractiviteSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type distances(distancesSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type placeslibres(placeslibresSEXP);
    rcpp_result_gen = Rcpp::wrap(one_distrib_continu(entrants, fuite, attractivite, distances, placeslibres));
    return rcpp_result_gen;
END_RCPP
}
// repartir_continu
std::vector<double> repartir_continu(const double actifs, const double fuite, const std::vector<double>& attractivite, const std::vector<double>& distances, std::vector<double>& placeslibres);
RcppExport SEXP _rmeaps_repartir_continu(SEXP actifsSEXP, SEXP fuiteSEXP, SEXP attractiviteSEXP, SEXP distancesSEXP, SEXP placeslibresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< const double >::type fuite(fuiteSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type attractivite(attractiviteSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type distances(distancesSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type placeslibres(placeslibresSEXP);
    rcpp_result_gen = Rcpp::wrap(repartir_continu(actifs, fuite, attractivite, distances, placeslibres));
    return rcpp_result_gen;
END_RCPP
}
// meaps_tension
List meaps_tension(IntegerMatrix rkdist, NumericVector emplois, NumericVector actifs, NumericMatrix modds, NumericVector f, IntegerMatrix shuf, std::string mode, Nullable<NumericVector> oddssubjectifs, int nthreads, bool progress, bool normalisation, double fuite_min, double seuil_newton, double seuil_dispo);
RcppExport SEXP _rmeaps_meaps_tension(SEXP rkdistSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP moddsSEXP, SEXP fSEXP, SEXP shufSEXP, SEXP modeSEXP, SEXP oddssubjectifsSEXP, SEXP nthreadsSEXP, SEXP progressSEXP, SEXP normalisationSEXP, SEXP fuite_minSEXP, SEXP seuil_newtonSEXP, SEXP seuil_dispoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type rkdist(rkdistSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type modds(moddsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type shuf(shufSEXP);
    Rcpp::traits::input_parameter< std::string >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type oddssubjectifs(oddssubjectifsSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< bool >::type normalisation(normalisationSEXP);
    Rcpp::traits::input_parameter< double >::type fuite_min(fuite_minSEXP);
    Rcpp::traits::input_parameter< double >::type seuil_newton(seuil_newtonSEXP);
    Rcpp::traits::input_parameter< double >::type seuil_dispo(seuil_dispoSEXP);
    rcpp_result_gen = Rcpp::wrap(meaps_tension(rkdist, emplois, actifs, modds, f, shuf, mode, oddssubjectifs, nthreads, progress, normalisation, fuite_min, seuil_newton, seuil_dispo));
    return rcpp_result_gen;
END_RCPP
}
// meaps_continu_test
NumericVector meaps_continu_test(IntegerVector j_dist, IntegerVector p_dist, NumericVector x_dist, NumericVector emplois, NumericVector actifs, NumericVector f, IntegerMatrix shuf, std::string attraction, double alpha, double beta, int nthreads, bool progress, bool normalisation, double fuite_min);
RcppExport SEXP _rmeaps_meaps_continu_test(SEXP j_distSEXP, SEXP p_distSEXP, SEXP x_distSEXP, SEXP emploisSEXP, SEXP actifsSEXP, SEXP fSEXP, SEXP shufSEXP, SEXP attractionSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP nthreadsSEXP, SEXP progressSEXP, SEXP normalisationSEXP, SEXP fuite_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type j_dist(j_distSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p_dist(p_distSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x_dist(x_distSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type emplois(emploisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type actifs(actifsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type shuf(shufSEXP);
    Rcpp::traits::input_parameter< std::string >::type attraction(attractionSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    Rcpp::traits::input_parameter< bool >::type normalisation(normalisationSEXP);
    Rcpp::traits::input_parameter< double >::type fuite_min(fuite_minSEXP);
    rcpp_result_gen = Rcpp::wrap(meaps_continu_test(j_dist, p_dist, x_dist, emplois, actifs, f, shuf, attraction, alpha, beta, nthreads, progress, normalisation, fuite_min));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_RankedRowSparseMatrix();

static const R_CallMethodDef CallEntries[] = {
    {"_rmeaps_chances_absorption", (DL_FUNC) &_rmeaps_chances_absorption, 4},
    {"_rmeaps_communaliser", (DL_FUNC) &_rmeaps_communaliser, 3},
    {"_rmeaps_meaps_oneshuf", (DL_FUNC) &_rmeaps_meaps_oneshuf, 11},
    {"_rmeaps_meaps_continu_cpp", (DL_FUNC) &_rmeaps_meaps_continu_cpp, 16},
    {"_rmeaps_meaps_multishuf", (DL_FUNC) &_rmeaps_meaps_multishuf, 13},
    {"_rmeaps_meaps_optim_cpp", (DL_FUNC) &_rmeaps_meaps_optim_cpp, 18},
    {"_rmeaps_one_distrib_continu", (DL_FUNC) &_rmeaps_one_distrib_continu, 5},
    {"_rmeaps_repartir_continu", (DL_FUNC) &_rmeaps_repartir_continu, 5},
    {"_rmeaps_meaps_tension", (DL_FUNC) &_rmeaps_meaps_tension, 14},
    {"_rmeaps_meaps_continu_test", (DL_FUNC) &_rmeaps_meaps_continu_test, 14},
    {"_rcpp_module_boot_RankedRowSparseMatrix", (DL_FUNC) &_rcpp_module_boot_RankedRowSparseMatrix, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_rmeaps(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
