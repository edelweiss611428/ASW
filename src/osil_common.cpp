#include "osil_common.h"

#include <cmath>

namespace osil {

std::vector<Index> clusterSizes(const Rcpp::IntegerVector& labels, int k) {
  std::vector<Index> size(k, 0);
  for (Index i = 0; i < labels.size(); ++i) {
    const int c = labels[i];
    if (c < 0 || c >= k) Rcpp::stop("cluster label out of range [0, k)");
    ++size[c];
  }
  return size;
}

void requirePartition(const std::vector<Index>& size) {
  for (std::size_t c = 0; c < size.size(); ++c) {
    if (size[c] == 0) Rcpp::stop("initial clustering contains an empty cluster");
  }
}

std::vector<double> clusterTotals(const Rcpp::IntegerVector& labels,
                                  const Rcpp::NumericVector& dist, Index n,
                                  int k) {
  std::vector<double> phi(static_cast<std::size_t>(n) * k, 0.0);
  for (Index j = 0; j + 1 < n; ++j) {
    Index pos = distIndex(n, j + 1, j);
    for (Index i = j + 1; i < n; ++i, ++pos) {
      const double d = dist[pos];
      phi[i * k + labels[j]] += d;
      phi[j * k + labels[i]] += d;
    }
  }
  return phi;
}

double aswFromLabels(const Rcpp::IntegerVector& labels,
                     const Rcpp::NumericVector& dist, Index n, int k) {
  const std::vector<Index> size = clusterSizes(labels, k);
  const std::vector<double> phi = clusterTotals(labels, dist, n, k);

  double total = 0.0;
  for (Index i = 0; i < n; ++i) {
    const int own = labels[i];
    if (size[own] < 2) continue;

    double b = kInf;
    for (int c = 0; c < k; ++c) {
      if (c == own || size[c] == 0) continue;
      b = std::min(b, phi[i * k + c] / size[c]);
    }
    if (!std::isfinite(b)) continue;  // no other non-empty cluster
    total += silhouette(phi[i * k + own] / (size[own] - 1), b);
  }
  return total / n;
}

}  // namespace osil
