#include "osil_core.h"

using namespace Rcpp;

// [[Rcpp::export(.effOSilCpp)]]
List effOSilCpp(NumericVector dist, IntegerVector iC, int N, int k) {
  osil::OSilFit fit = osil::effOSilCore(dist, iC, N, k);
  return List::create(_["Clustering"] = fit.labels + 1L, _["ASW"] = fit.asw,
                      _["nIter"] = fit.iterations);
}

// [[Rcpp::export(.OSilCpp)]]
List OSilCpp(NumericVector dist, IntegerVector iC, int N, int k) {
  osil::OSilNaiveFit fit = osil::osilCore(dist, iC, N, k);
  return List::create(_["Clustering"] = fit.labels + 1L, _["ASW"] = fit.asw,
                      _["nIter"] = fit.iterations);
}
