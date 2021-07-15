// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// MoransI
double MoransI(NumericMatrix Mat1, int r1);
RcppExport SEXP _ICvectorfields_MoransI(SEXP Mat1SEXP, SEXP r1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Mat1(Mat1SEXP);
    Rcpp::traits::input_parameter< int >::type r1(r1SEXP);
    rcpp_result_gen = Rcpp::wrap(MoransI(Mat1, r1));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ICvectorfields_MoransI", (DL_FUNC) &_ICvectorfields_MoransI, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ICvectorfields(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
