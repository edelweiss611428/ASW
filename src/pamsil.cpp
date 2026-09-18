#include "osil_common.h"

using namespace Rcpp;
using osil::Index;
using osil::kInf;

namespace {

double medoidDist(const NumericVector& dist, Index n, Index i, Index m) {
  return i == m ? 0.0 : dist[osil::distIndex(n, i, m)];
}

// Distance to the own medoid and to the nearest other medoid, with its label.
void nearestTwo(const NumericVector& dist, const IntegerVector& medoids,
                const IntegerVector& labels, IntegerVector& second,
                NumericVector& d1, NumericVector& d2, Index n, int k) {
  for (Index i = 0; i < n; ++i) {
    const int own = labels[i];
    double best = kInf;
    int bestLab = own;

    for (int c = 0; c < k; ++c) {
      if (c == own) continue;
      const double d = medoidDist(dist, n, i, medoids[c]);
      if (d < best) {
        best = d;
        bestLab = c;
      }
    }

    second[i] = bestLab;
    d2[i] = best;
    d1[i] = medoidDist(dist, n, i, medoids[own]);
  }
}

// Clustering obtained by making point newMedoid the medoid of cluster target,
// in O(n): only one medoid changes, so each point compares its distance to the
// new medoid against its current best (or second best, if its own medoid went).
IntegerVector assignAfterSwap(const NumericVector& dist,
                              const IntegerVector& medoids,
                              const IntegerVector& labels,
                              const IntegerVector& second,
                              const NumericVector& d1, const NumericVector& d2,
                              Index n, Index newMedoid, int target) {
  IntegerVector out(n);

  for (Index i = 0; i < n; ++i) {
    const int own = labels[i];

    if (i == newMedoid) {
      out[i] = target;
    } else if (own != target && i == medoids[own]) {
      out[i] = own;
    } else {
      const double td = medoidDist(dist, n, i, newMedoid);
      if (own != target) {
        out[i] = td < d1[i] ? target : own;
      } else {  // the medoid of this point's own cluster has been replaced
        out[i] = td < d2[i] ? target : second[i];
      }
    }
  }
  return out;
}

}  // namespace

// [[Rcpp::export(.PAMSilCpp)]]
List PAMSilCpp(NumericVector dist, IntegerVector C, IntegerVector medoids,
               int N, int k) {
  IntegerVector labels = clone(C);
  IntegerVector med = clone(medoids);

  IntegerVector second(N);
  NumericVector d1(N), d2(N);
  nearestTwo(dist, med, labels, second, d1, d2, N, k);

  double best = osil::aswFromLabels(labels, dist, N, k);
  int iterations = 0;
  bool moved = true;

  while (moved) {
    ++iterations;
    moved = false;
    Index bestPoint = 0;
    int bestCluster = 0;

    for (Index i = 0; i < N; ++i) {
      if (i == med[labels[i]]) continue;
      for (int j = 0; j < k; ++j) {
        const IntegerVector trial =
            assignAfterSwap(dist, med, labels, second, d1, d2, N, i, j);
        const double candidate = osil::aswFromLabels(trial, dist, N, k);
        if (candidate > best) {
          best = candidate;
          bestPoint = i;
          bestCluster = j;
          moved = true;
        }
      }
      checkUserInterrupt();
    }

    if (moved) {
      labels = assignAfterSwap(dist, med, labels, second, d1, d2, N, bestPoint,
                               bestCluster);
      med[bestCluster] = bestPoint;
      nearestTwo(dist, med, labels, second, d1, d2, N, k);
    }
  }

  return List::create(_["Clustering"] = labels + 1L, _["medoids"] = med + 1L,
                      _["ASW"] = best, _["nIter"] = iterations);
}
