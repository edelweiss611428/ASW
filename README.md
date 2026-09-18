# ASW

Clustering by optimising the Average Silhouette Width.

The ASW is usually used to *validate* a clustering. It can also be optimised
directly, which Batool & Hennig (2021) showed produces good clusterings — but
their OSil algorithm costs `O(k N^3)` per iteration, which puts it out of reach
well before most people expect. 

This package implements `effOSil`, which returns the *same* clustering as OSil
with an `O(N)` reduction in cost, and `scalOSil`, which does the same for FOSil.

| Function | Use when |
|---|---|
| `effOSil()` | the distance matrix fits in memory |
| `scalOSil()` | it does not, or N is large enough that `effOSil()` is still slow |
| `PAMSil()` | you want the medoid-based objective of Van der Laan et al. (2003) |
| `Init()` | you want the best of several standard clusterings by ASW |
| `asw()`, `Silhouette()` | you have a clustering and want to score it |

```r
library(ASW)

dx  = dist(scale(faithful))
fit = effOSil(dx, K = 2:12)

fit$k                  # selected number of clusters
fit$best_asw           # ASW attained
plot(faithful, col = fit$best_clustering)
```

Both `effOSil()` and `scalOSil()` take a `variant` argument selecting the
original algorithm instead, for timing comparisons:

```r
system.time(effOSil(dx, K = 5, variant = "efficient"))
system.time(effOSil(dx, K = 5, variant = "original"))
```

## Installation

```r
# install.packages("remotes")
remotes::install_github("edelweiss611428/ASW")
```

## Notes

Observations in singleton clusters are given a silhouette width of 0, following
Rousseeuw (1987), as are coincident observations for which both the within- and
between-cluster mean distances are zero.

`effOSil()` and OSil agree exactly on data without duplicated observations. When
the data contains duplicates the ASW has exact ties, and the two may resolve a
tied reassignment differently and converge to different local optima of equal or
near-equal ASW. The same applies to `scalOSil()` and FOSil.

## References

Batool, F. and Hennig, C. (2021). Clustering with the average silhouette width.
*Computational Statistics & Data Analysis*, 158, 107190.

Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the interpretation and
validation of cluster analysis. *Journal of Computational and Applied
Mathematics*, 20, 53–65.

Van der Laan, M., Pollard, K. and Bryan, J. (2003). A new partitioning around
medoids algorithm. *Journal of Statistical Computation and Simulation*, 73(8),
575–584.
