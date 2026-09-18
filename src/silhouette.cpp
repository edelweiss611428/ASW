#include "osil_common.h"

#include <cmath>

using namespace Rcpp;
using osil::Index;

// [[Rcpp::export(.ASWCpp)]]
double ASWCpp(IntegerVector C, NumericVector dist, int N, int k) {
  return osil::aswFromLabels(C, dist, N, k);
}

// [[Rcpp::export(.SWCpp)]]
List SWCpp(IntegerVector C, NumericVector dist, int N, int k) {
  const std::vector<Index> size = osil::clusterSizes(C, k);
  const std::vector<double> phi = osil::clusterTotals(C, dist, N, k);

  NumericVector width(N);
  IntegerVector neighbour(N, NA_INTEGER);

  for (Index i = 0; i < N; ++i) {
    const int own = C[i];

    double b = osil::kInf;
    int lb = NA_INTEGER;
    for (int c = 0; c < k; ++c) {
      if (c == own || size[c] == 0) continue;
      const double v = phi[i * k + c] / size[c];
      if (v < b) {
        b = v;
        lb = c + 1;
      }
    }

    if (!std::isfinite(b)) continue;  // no other non-empty cluster
    neighbour[i] = lb;                // well defined even for a singleton
    if (size[own] < 2) continue;      // but a singleton has width 0

    width[i] = osil::silhouette(phi[i * k + own] / (size[own] - 1), b);
  }

  return List::create(_["neighbor"] = neighbour, _["sil_width"] = width);
}

// [[Rcpp::export(.subDistCpp)]]
NumericVector subDistCpp(NumericVector dist, IntegerVector idx, bool diag,
                         bool upper, int N, int n) {
  NumericVector sub(static_cast<Index>(n) * (n - 1) / 2);

  Index l = 0;
  for (Index j = 0; j + 1 < n; ++j) {
    for (Index i = j + 1; i < n; ++i) {
      sub[l++] = dist[osil::distIndex(N, idx[i], idx[j])];
    }
  }

  sub.attr("Size") = n;
  sub.attr("Diag") = diag;
  sub.attr("Upper") = upper;
  sub.attr("class") = "dist";
  return sub;
}
