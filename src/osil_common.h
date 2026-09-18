#ifndef OSIL_COMMON_H
#define OSIL_COMMON_H

#include <Rcpp.h>

#include <algorithm>
#include <limits>
#include <vector>

namespace osil {

using Index = R_xlen_t;

constexpr double kInf = std::numeric_limits<double>::infinity();

// Position of d(i, j) in a "dist" object over n points. 0-based, requires i != j.
inline Index distIndex(Index n, Index i, Index j) noexcept {
  if (i < j) std::swap(i, j);
  return (2 * n - 1 - j) * j / 2 + (i - j) - 1;
}

inline double silhouette(double a, double b) noexcept {
  const double m = std::max(a, b);
  return m > 0.0 ? (b - a) / m : 0.0;
}

std::vector<Index> clusterSizes(const Rcpp::IntegerVector& labels, int k);

void requirePartition(const std::vector<Index>& size);

// phi[i * k + c] = sum of distances from point i to the members of cluster c.
std::vector<double> clusterTotals(const Rcpp::IntegerVector& labels,
                                  const Rcpp::NumericVector& dist, Index n,
                                  int k);

double aswFromLabels(const Rcpp::IntegerVector& labels,
                     const Rcpp::NumericVector& dist, Index n, int k);

}  // namespace osil

#endif  // OSIL_COMMON_H
