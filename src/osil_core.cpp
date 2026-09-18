#include "osil_core.h"

using namespace Rcpp;

namespace osil {

void NeighbourCache::refresh(const IntegerVector& labels,
                             const std::vector<double>& phi,
                             const std::vector<Index>& size, Index n, int k) {
  for (Index i = 0; i < n; ++i) {
    const int own = labels[i];
    double b1 = kInf, b2 = kInf, b3 = kInf;
    int l1 = -1, l2 = -1, l3 = -1;

    for (int c = 0; c < k; ++c) {
      if (c == own || size[c] == 0) continue;
      const double v = phi[i * k + c] / size[c];
      if (v < b1) {
        b3 = b2; l3 = l2;
        b2 = b1; l2 = l1;
        b1 = v;  l1 = c;
      } else if (v < b2) {
        b3 = b2; l3 = l2;
        b2 = v;  l2 = c;
      } else if (v < b3) {
        b3 = v;  l3 = c;
      }
    }

    Neighbours& r = row[i];
    r.b = b1;  r.s = b2;  r.h = b3;
    r.lb = l1; r.ls = l2; r.lh = l3;
    r.a = size[own] > 1 ? phi[i * k + own] / (size[own] - 1) : kInf;
  }
}

double NeighbourCache::asw(const IntegerVector& labels,
                           const std::vector<Index>& size, Index n) const {
  double total = 0.0;
  for (Index i = 0; i < n; ++i) {
    // singletons, and points with no other non-empty cluster, contribute 0
    if (size[labels[i]] < 2 || !std::isfinite(row[i].b)) continue;
    total += silhouette(row[i].a, row[i].b);
  }
  return total / n;
}

namespace {

// Distances from point i to every other point, gathered into drow. Half of a
// "dist" row is strided, so this is done once per i and reused for all k
// candidate clusters rather than re-gathered inside the j loop.
void gatherRow(const NumericVector& dist, Index n, Index i, double* drow) {
  drow[i] = 0.0;
  for (Index l = 0; l < i; ++l) drow[l] = dist[distIndex(n, i, l)];
  if (i + 1 < n) {
    Index pos = distIndex(n, i + 1, i);
    for (Index l = i + 1; l < n; ++l) drow[l] = dist[pos++];
  }
}

// ASW of every move of point i, in one O(n) pass. cand[j] receives the ASW of
// moving i into cluster j (cand[from] is left at 0 and must be ignored).
// Evaluating all k targets together is what makes this cheap: everything that
// does not depend on j -- the loads of labels[l], drow[l] and nc.row[l], the
// own == from test, and the divisions by size[from] and size[own] -- is done
// once per l rather than k - 1 times.
// Requires size[labels[i]] >= 2 and drow filled by gatherRow.
void sweepPoint(const double* drow, const IntegerVector& labels,
                const std::vector<double>& phi, const std::vector<Index>& size,
                const NeighbourCache& nc, Index n, int k, Index i, int from,
                double* cand) {
  for (int j = 0; j < k; ++j) cand[j] = 0.0;
  const Neighbours& ri = nc.row[i];

  for (Index l = 0; l < n; ++l) {
    const double* pl = &phi[l * k];

    if (l == i) {  // the moving point: its old cluster is now at distance ri.a
      for (int j = 0; j < k; ++j) {
        if (j == from) continue;
        cand[j] += silhouette(pl[j] / size[j],
                              std::min(ri.a, ri.lb == j ? ri.s : ri.b));
      }
      continue;
    }

    const int own = labels[l];
    const double td = drow[l];
    const Neighbours& r = nc.row[l];

    if (own == from) {  // l stays behind in the shrinking cluster
      if (size[from] == 2) continue;  // l is left alone, silhouette 0
      const double a = (pl[from] - td) / (size[from] - 2);
      for (int j = 0; j < k; ++j) {
        if (j == from) continue;
        const double grown = (pl[j] + td) / (size[j] + 1);
        cand[j] += silhouette(a, std::min(grown, r.lb == j ? r.s : r.b));
      }
    } else {
      // A singleton stays at silhouette 0 unless it is the cluster receiving i,
      // so this is skipped per j, not for l as a whole.
      const bool singleton = (size[own] == 1);
      const double shrunk = (pl[from] - td) / (size[from] - 1);

      // Nearest cluster other than l's own, from, and j. As j varies this
      // takes one value for every j but a single exceptional index.
      double restDef, restExc;
      int excJ;
      if (r.lb != from) {
        restDef = r.b; excJ = r.lb; restExc = (r.ls == from) ? r.h : r.s;
      } else {
        restDef = r.s; excJ = r.ls; restExc = r.h;
      }
      // j == own: l is in the receiving cluster. Independent of j given own.
      const double aOwn = (pl[own] + td) / size[own];
      const double bOwn = std::min(shrunk, r.lb == from ? r.s : r.b);

      for (int j = 0; j < k; ++j) {
        if (j == from) continue;
        if (j == own) {
          cand[j] += silhouette(aOwn, bOwn);
          continue;
        }
        if (singleton) continue;
        const double grown = (pl[j] + td) / (size[j] + 1);
        const double rest = (j == excJ) ? restExc : restDef;
        cand[j] += silhouette(r.a, std::min(std::min(shrunk, grown), rest));
      }
    }
  }
  for (int j = 0; j < k; ++j) cand[j] /= n;
}

void applyMove(const NumericVector& dist, std::vector<double>& phi, Index n,
               int k, Index i, int from, int to) {
  for (Index l = 0; l < n; ++l) {
    if (l == i) continue;
    const double td = dist[distIndex(n, i, l)];
    phi[l * k + from] -= td;
    phi[l * k + to] += td;
  }
}

}  // namespace

OSilFit effOSilCore(const NumericVector& dist, const IntegerVector& init,
                    Index n, int k) {
  IntegerVector labels = clone(init);
  std::vector<Index> size = clusterSizes(labels, k);
  requirePartition(size);
  std::vector<double> phi = clusterTotals(labels, dist, n, k);

  NeighbourCache cache(n);
  cache.refresh(labels, phi, size, n, k);
  double best = cache.asw(labels, size, n);

  std::vector<double> drow(n), cand(k);
  int iterations = 0;
  bool moved = true;

  while (moved) {
    ++iterations;
    moved = false;
    Index bestPoint = 0;
    int bestCluster = 0;

    for (Index i = 0; i < n; ++i) {
      const int from = labels[i];
      if (size[from] == 1) continue;  // never empty a cluster
      gatherRow(dist, n, i, drow.data());
      sweepPoint(drow.data(), labels, phi, size, cache, n, k, i, from, cand.data());
      for (int j = 0; j < k; ++j) {
        if (j == from) continue;
        if (cand[j] > best) {
          best = cand[j];
          bestPoint = i;
          bestCluster = j;
          moved = true;
        }
      }
    }

    if (moved) {
      const int from = labels[bestPoint];
      applyMove(dist, phi, n, k, bestPoint, from, bestCluster);
      labels[bestPoint] = bestCluster;
      --size[from];
      ++size[bestCluster];
      cache.refresh(labels, phi, size, n, k);
    }
    checkUserInterrupt();
  }

  // phi is maintained in place across the run, so it carries O(q * eps * max|phi|)
  // of accumulated rounding. Recomputing the final ASW from the labels costs
  // O(n^2) once -- negligible against O(q k n^2) -- and makes the reported value
  // exact for the reported clustering.
  best = aswFromLabels(labels, dist, n, k);

  return OSilFit{labels, best, iterations, std::move(phi), std::move(size),
                 std::move(cache)};
}

OSilNaiveFit osilCore(const NumericVector& dist, const IntegerVector& init,
                      Index n, int k) {
  IntegerVector labels = clone(init);
  std::vector<Index> size = clusterSizes(labels, k);
  requirePartition(size);

  double best = aswFromLabels(labels, dist, n, k);
  IntegerVector trial = clone(labels);

  int iterations = 0;
  bool moved = true;

  while (moved) {
    ++iterations;
    moved = false;
    Index bestPoint = 0;
    int bestCluster = 0;

    for (Index i = 0; i < n; ++i) {
      const int from = labels[i];
      if (size[from] == 1) continue;
      for (int j = 0; j < k; ++j) {
        if (j == from) continue;
        trial[i] = j;
        const double candidate = aswFromLabels(trial, dist, n, k);
        if (candidate > best) {
          best = candidate;
          bestPoint = i;
          bestCluster = j;
          moved = true;
        }
      }
      trial[i] = from;
    }

    if (moved) {
      --size[labels[bestPoint]];
      ++size[bestCluster];
      labels[bestPoint] = bestCluster;
      trial[bestPoint] = bestCluster;
    }
    checkUserInterrupt();
  }

  return OSilNaiveFit{labels, best, iterations};
}

}  // namespace osil
