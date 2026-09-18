#include "osil_core.h"

using namespace Rcpp;
using osil::Index;
using osil::kInf;

// [[Rcpp::export(.scalOSil_PC)]]
List scalOSil_PC(NumericVector dist, IntegerVector iC, int N, int k) {
  osil::OSilFit fit = osil::effOSilCore(dist, iC, N, k);

  NumericVector ai(N), bi(N), si(N);
  IntegerVector lbi(N);
  for (Index i = 0; i < N; ++i) {
    const osil::Neighbours& r = fit.cache.row[i];
    ai[i] = r.a; bi[i] = r.b; si[i] = r.s; lbi[i] = r.lb;
  }

  return List::create(
      _["PC"] = fit.labels, _["ASW"] = fit.asw,
      _["PC_UCT"] = NumericVector(fit.phi.begin(), fit.phi.end()),
      _["Nj"] = IntegerVector(fit.size.begin(), fit.size.end()),
      _["PC_ai"] = ai, _["PC_bi"] = bi, _["PC_si"] = si, _["PC_lbi"] = lbi);
}

// Classify each unassigned point into the cluster of the partial clustering
// that maximises the ASW of PC + {y}, each point independently. The constant
// 1 / (n1 + 1) of the ASW is dropped: it does not affect the argmax.
// [[Rcpp::export(.scalOSil_C)]]
IntegerVector scalOSil_C(NumericVector dist, int k, List PC_result,
                         IntegerVector idxPC, IntegerVector idxC, int n1,
                         int n2, int N) {
  const IntegerVector pc = PC_result["PC"];
  const NumericVector phi = PC_result["PC_UCT"];
  const IntegerVector size = PC_result["Nj"];
  const NumericVector ca = PC_result["PC_ai"];
  const NumericVector cb = PC_result["PC_bi"];
  const NumericVector cs = PC_result["PC_si"];
  const IntegerVector clb = PC_result["PC_lbi"];

  IntegerVector out(N);
  std::copy(pc.begin(), pc.end(), out.begin());

  std::vector<double> td(n1), toCluster(k);

  for (int i = 0; i < n2; ++i) {
    const Index y = idxC[i];
    std::fill(toCluster.begin(), toCluster.end(), 0.0);

    for (int j = 0; j < n1; ++j) {
      const double d = dist[osil::distIndex(N, y, idxPC[j])];
      td[j] = d;
      toCluster[pc[j]] += d;
    }

    double best = -kInf;
    int bestCluster = 0;

    for (int m = 0; m < k; ++m) {
      if (size[m] == 0) continue;
      double total = 0.0;

      for (int j = 0; j < n1; ++j) {
        const int own = pc[j];
        double a, b;
        if (own == m) {  // y joins this point's own cluster; b is unchanged
          a = (phi[j * k + m] + td[j]) / size[m];
          b = cb[j];
        } else {
          if (size[own] == 1) continue;
          a = ca[j];
          const double grown = (phi[j * k + m] + td[j]) / (size[m] + 1);
          b = std::min(grown, clb[j] == m ? cs[j] : cb[j]);
        }
        total += osil::silhouette(a, b);
      }

      double ay = toCluster[m] / size[m], by = kInf;
      for (int c = 0; c < k; ++c) {
        if (c == m || size[c] == 0) continue;
        by = std::min(by, toCluster[c] / size[c]);
      }
      total += osil::silhouette(ay, by);

      if (total > best) {
        best = total;
        bestCluster = m;
      }
    }

    out[n1 + i] = bestCluster;
    checkUserInterrupt();
  }
  return out;
}
