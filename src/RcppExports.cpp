// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// grad_logistic
arma::mat grad_logistic(arma::vec& x, List opts);
RcppExport SEXP _apg_grad_logistic(SEXP xSEXP, SEXP optsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type opts(optsSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_logistic(x, opts));
    return rcpp_result_gen;
END_RCPP
}
// grad_ranking_logistic
NumericVector grad_ranking_logistic(NumericVector x, List opts);
RcppExport SEXP _apg_grad_ranking_logistic(SEXP xSEXP, SEXP optsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type opts(optsSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_ranking_logistic(x, opts));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_apg_grad_logistic", (DL_FUNC) &_apg_grad_logistic, 2},
    {"_apg_grad_ranking_logistic", (DL_FUNC) &_apg_grad_ranking_logistic, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_apg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
