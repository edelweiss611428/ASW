// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ASWCpp
double ASWCpp(IntegerVector C, NumericVector dist, int N, int k);
RcppExport SEXP _EfficientOASW_ASWCpp(SEXP CSEXP, SEXP distSEXP, SEXP NSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(ASWCpp(C, dist, N, k));
    return rcpp_result_gen;
END_RCPP
}
// subDistCpp
NumericVector subDistCpp(NumericVector dist, IntegerVector idx, bool diag, bool upper, int N, int n);
RcppExport SEXP _EfficientOASW_subDistCpp(SEXP distSEXP, SEXP idxSEXP, SEXP diagSEXP, SEXP upperSEXP, SEXP NSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< bool >::type diag(diagSEXP);
    Rcpp::traits::input_parameter< bool >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(subDistCpp(dist, idx, diag, upper, N, n));
    return rcpp_result_gen;
END_RCPP
}
// SWCpp
List SWCpp(IntegerVector C, NumericVector dist, int N, int k);
RcppExport SEXP _EfficientOASW_SWCpp(SEXP CSEXP, SEXP distSEXP, SEXP NSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(SWCpp(C, dist, N, k));
    return rcpp_result_gen;
END_RCPP
}
// effOSilCpp
List effOSilCpp(NumericVector dist, IntegerVector iC, int N, int k);
RcppExport SEXP _EfficientOASW_effOSilCpp(SEXP distSEXP, SEXP iCSEXP, SEXP NSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type iC(iCSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(effOSilCpp(dist, iC, N, k));
    return rcpp_result_gen;
END_RCPP
}
// OSilCpp
List OSilCpp(NumericVector dist, IntegerVector iC, int N, int k);
RcppExport SEXP _EfficientOASW_OSilCpp(SEXP distSEXP, SEXP iCSEXP, SEXP NSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type iC(iCSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(OSilCpp(dist, iC, N, k));
    return rcpp_result_gen;
END_RCPP
}
// scalOSil_PC_Step
List scalOSil_PC_Step(NumericVector dist, IntegerVector iC, int N, int k);
RcppExport SEXP _EfficientOASW_scalOSil_PC_Step(SEXP distSEXP, SEXP iCSEXP, SEXP NSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type iC(iCSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(scalOSil_PC_Step(dist, iC, N, k));
    return rcpp_result_gen;
END_RCPP
}
// scalOSil_C_Step
IntegerVector scalOSil_C_Step(NumericVector dist, int k, List PC_result, IntegerVector idxPC, IntegerVector idxC, int n1, int n2, int N);
RcppExport SEXP _EfficientOASW_scalOSil_C_Step(SEXP distSEXP, SEXP kSEXP, SEXP PC_resultSEXP, SEXP idxPCSEXP, SEXP idxCSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< List >::type PC_result(PC_resultSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idxPC(idxPCSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idxC(idxCSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(scalOSil_C_Step(dist, k, PC_result, idxPC, idxC, n1, n2, N));
    return rcpp_result_gen;
END_RCPP
}
// FOSil_PC_Step
List FOSil_PC_Step(NumericVector dist, IntegerVector iC, int N, int k);
RcppExport SEXP _EfficientOASW_FOSil_PC_Step(SEXP distSEXP, SEXP iCSEXP, SEXP NSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type iC(iCSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(FOSil_PC_Step(dist, iC, N, k));
    return rcpp_result_gen;
END_RCPP
}
// FOSil_C_Step
IntegerVector FOSil_C_Step(NumericVector dist, NumericVector distPC, int k, List PC_result, IntegerVector idxPC, IntegerVector idxC, int n1, int n2, int N);
RcppExport SEXP _EfficientOASW_FOSil_C_Step(SEXP distSEXP, SEXP distPCSEXP, SEXP kSEXP, SEXP PC_resultSEXP, SEXP idxPCSEXP, SEXP idxCSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type distPC(distPCSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< List >::type PC_result(PC_resultSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idxPC(idxPCSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idxC(idxCSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(FOSil_C_Step(dist, distPC, k, PC_result, idxPC, idxC, n1, n2, N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EfficientOASW_ASWCpp", (DL_FUNC) &_EfficientOASW_ASWCpp, 4},
    {"_EfficientOASW_subDistCpp", (DL_FUNC) &_EfficientOASW_subDistCpp, 6},
    {"_EfficientOASW_SWCpp", (DL_FUNC) &_EfficientOASW_SWCpp, 4},
    {"_EfficientOASW_effOSilCpp", (DL_FUNC) &_EfficientOASW_effOSilCpp, 4},
    {"_EfficientOASW_OSilCpp", (DL_FUNC) &_EfficientOASW_OSilCpp, 4},
    {"_EfficientOASW_scalOSil_PC_Step", (DL_FUNC) &_EfficientOASW_scalOSil_PC_Step, 4},
    {"_EfficientOASW_scalOSil_C_Step", (DL_FUNC) &_EfficientOASW_scalOSil_C_Step, 8},
    {"_EfficientOASW_FOSil_PC_Step", (DL_FUNC) &_EfficientOASW_FOSil_PC_Step, 4},
    {"_EfficientOASW_FOSil_C_Step", (DL_FUNC) &_EfficientOASW_FOSil_C_Step, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_EfficientOASW(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
