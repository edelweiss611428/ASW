#ifndef OSIL_CORE_H
#define OSIL_CORE_H

#include "osil_common.h"

#include <cmath>

namespace osil {

// Mean distance from a point to its own cluster (a) and to the three nearest
// other clusters (b <= s <= h), with their labels. Moving one point changes at
// most two clusters, so three candidates always leave one untouched -- that is
// exactly why b, s, h suffice and no more. Absent entries are +Inf / -1, which
// removes the k = 2 and k = 3 special cases.
// Kept as one struct per point rather than seven parallel arrays: the update
// loop reads all seven fields for the same point, so this is one cache line
// instead of seven streams.
struct Neighbours {
  double a = kInf, b = kInf, s = kInf, h = kInf;
  int lb = -1, ls = -1, lh = -1;
};

struct NeighbourCache {
  std::vector<Neighbours> row;

  explicit NeighbourCache(Index n) : row(n) {}

  void refresh(const Rcpp::IntegerVector& labels, const std::vector<double>& phi,
               const std::vector<Index>& size, Index n, int k);

  double asw(const Rcpp::IntegerVector& labels, const std::vector<Index>& size,
             Index n) const;
};

struct OSilFit {
  Rcpp::IntegerVector labels;  // 0-based
  double asw;
  int iterations;
  std::vector<double> phi;
  std::vector<Index> size;
  NeighbourCache cache;
};

struct OSilNaiveFit {
  Rcpp::IntegerVector labels;  // 0-based
  double asw;
  int iterations;
};

// Steepest-ascent local search: one sweep evaluates every (point, cluster)
// move and applies the single best improving one. effOSilCore scores a move in
// O(n) from the cache; osilCore rebuilds the ASW from scratch in O(n^2).
OSilFit effOSilCore(const Rcpp::NumericVector& dist,
                    const Rcpp::IntegerVector& init, Index n, int k);

OSilNaiveFit osilCore(const Rcpp::NumericVector& dist,
                      const Rcpp::IntegerVector& init, Index n, int k);

}  // namespace osil

#endif  // OSIL_CORE_H
