#include "osil_core.h"

using namespace Rcpp;
using osil::Index;
using osil::kInf;

namespace {

// Grow a dist object over n points into one over n + 1; the new row is zero.
NumericVector padDist(const NumericVector& d, Index n) {
  NumericVector out((n + 1) * n / 2);
  Index l = 0;
  for (Index j = 0; j < n; ++j) {
    for (Index i = j + 1; i <= n; ++i) {
      out[l++] = (i == n) ? 0.0 : d[osil::distIndex(n, i, j)];
    }
  }
  return out;
}

}  // namespace

// [[Rcpp::export(.FOSil_PC)]]
List FOSil_PC(NumericVector dist, IntegerVector iC, int N, int k) {
  osil::OSilNaiveFit fit = osil::osilCore(dist, iC, N, k);
  const std::vector<Index> size = osil::clusterSizes(fit.labels, k);
  return List::create(_["PC"] = fit.labels, _["ASW"] = fit.asw,
                      _["Nj"] = IntegerVector(size.begin(), size.end()));
}

// [[Rcpp::export(.FOSil_C)]]
IntegerVector FOSil_C(NumericVector dist, NumericVector distPC, int k,
                      List PC_result, IntegerVector idxPC, IntegerVector idxC,
                      int n1, int n2, int N) {
  const IntegerVector pc = PC_result["PC"];

  IntegerVector out(N);
  std::copy(pc.begin(), pc.end(), out.begin());

  NumericVector padded = padDist(distPC, n1);
  IntegerVector trial(n1 + 1);
  std::copy(pc.begin(), pc.end(), trial.begin());

  for (int i = 0; i < n2; ++i) {
    for (int j = 0; j < n1; ++j) {
      padded[osil::distIndex(n1 + 1, n1, j)] =
          dist[osil::distIndex(N, idxPC[j], idxC[i])];
    }

    double best = -kInf;
    int bestCluster = 0;
    for (int m = 0; m < k; ++m) {
      trial[n1] = m;
      const double candidate = osil::aswFromLabels(trial, padded, n1 + 1, k);
      if (candidate > best) {
        best = candidate;
        bestCluster = m;
      }
    }

    out[n1 + i] = bestCluster;
    checkUserInterrupt();
  }
  return out;
}
