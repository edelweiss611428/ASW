#' @name scalOSil
#' @title The Scalable Optimum Silhouette algorithm
#'
#' @description Clusters a distance matrix by maximising the Average Silhouette
#' Width (ASW), using the Scalable Optimum Silhouette (scalOSil) algorithm.
#' Intended for data sets too large for [effOSil].
#'
#' @param dx A `dist` object, as returned by [stats::dist()].
#' @param K An integer vector specifying the numbers of clusters to consider.
#' Defaults to `2:12`.
#' @param n The subsample size. Defaults to `ceiling(0.1 * N)`, where `N` is the
#' number of observations.
#' @param ns The number of subsamples drawn per instance. Defaults to `10`.
#' @param rep The number of scalOSil instances. Defaults to `1`.
#' @param initMethod A character vector of initialisation methods. Defaults to
#' `"average"`. See [Init].
#' @param variant Either `"scalable"` (the default), which uses scalOSil, or
#' `"original"`, which uses the original FOSil algorithm. Both classify
#' observations by the same rule; `"original"` is provided for comparison and is
#' far slower.
#'
#' @return An object of class `"ASW"`, a list with components:
#' \describe{
#'   \item{best_clustering}{An integer vector giving the clustering that
#'     achieves the highest ASW.}
#'   \item{best_asw}{The highest ASW value.}
#'   \item{k}{The estimated number of clusters.}
#'   \item{clusterings}{An integer matrix with one column per element of `K`.}
#'   \item{asw}{A numeric vector of ASW values, one per element of `K`.}
#'   \item{n}{The subsample size used.}
#'   \item{method}{The name of the algorithm used.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details scalOSil follows the two-step structure of FOSil (Batool & Hennig,
#' 2021): a partial clustering step, which runs [effOSil] on a subsample of size
#' `n`, and a classification step, which assigns each remaining observation to
#' the cluster maximising the ASW of the partial clustering together with that
#' observation. Observations are classified independently of one another, so the
#' order of classification does not matter.
#'
#' The classification step returns exactly the same assignments as FOSil, but
#' evaluates each candidate in \eqn{O(n)} rather than \eqn{O(n^2)} time, giving
#' an \eqn{O(n)} reduction in its runtime.
#'
#' `ns` subsamples are drawn per instance and the one achieving the highest ASW
#' in the partial clustering step is carried forward. `rep` instances are run
#' and the one achieving the highest ASW overall is returned. This function uses
#' the random number generator; call [set.seed()] beforehand for reproducible
#' results.
#'
#' @examples
#' set.seed(59)
#' x = scale(faithful)
#' dx = dist(x)
#' fit = scalOSil(dx = dx, K = 2:8, n = ceiling(0.25 * nrow(x)), ns = 10, rep = 1)
#'
#' oldpar = par(mfrow = c(1, 2))
#' plot(faithful, col = fit$best_clustering, pch = fit$best_clustering)
#' plot(2:8, fit$asw, type = "l", xlab = "k", ylab = "ASW")
#' par(oldpar)
#'
#' @references
#' Batool, F. and Hennig, C. (2021). Clustering with the average silhouette
#' width. \emph{Computational Statistics & Data Analysis}, 158, 107190.
#' \doi{10.1016/j.csda.2021.107190}
#'
#' @seealso [effOSil], [PAMSil], [Init], [asw], [Silhouette].
#'
#' @author Minh Long Nguyen \email{edelweiss611428@@gmail.com}
#' @export

scalOSil = function(dx, K = 2:12, n = NULL, ns = 10, rep = 1,
                    initMethod = "average", variant = c("scalable", "original")){

  N = .checkDist(dx)

  if(is.null(n)){
    n = ceiling(0.1 * N)
  }

  n = .checkCount(n, "n", lower = 2L)

  if(n > N){
    stop("`n` cannot exceed the number of observations.", call. = FALSE)
  }

  ns = .checkCount(ns, "ns")
  rep = .checkCount(rep, "rep")
  K = .checkK(K, n, "the subsample size")
  initMethod = .checkInitMethod(initMethod)
  variant = match.arg(variant)

  if(variant == "scalable"){
    pcStep = function(dsub, init, k) .scalOSil_PC(dsub, init, n, k)
    cStep = function(dsub, k, pc, idxPC, idxC) .scalOSil_C(dx, k, pc, idxPC, idxC, n, N - n, N)
  } else {
    pcStep = function(dsub, init, k) .FOSil_PC(dsub, init, n, k)
    cStep = function(dsub, k, pc, idxPC, idxC) .FOSil_C(dx, dsub, k, pc, idxPC, idxC, n, N - n, N)
  }

  nK = length(K)
  clusterings = matrix(NA_integer_, nrow = N, ncol = nK, dimnames = list(NULL, K))
  asw = setNames(numeric(nK), K)

  for(i in seq_len(nK)){

    bestASW = -Inf
    bestClustering = NULL

    for(j in seq_len(rep)){

      bestPC = -Inf
      keep = NULL

      for(l in seq_len(ns)){

        idx = sample.int(N)
        idxPC = idx[seq_len(n)]
        dsub = .subDistCpp(dx, idxPC - 1L, FALSE, FALSE, N, n)
        init = Init(dsub, K[i], initMethod)$clustering - 1L
        pc = pcStep(dsub, init, K[i])

        if(pc$ASW > bestPC){
          bestPC = pc$ASW
          keep = list(idx = idx, idxPC = idxPC, idxC = idx[(n + 1L):N], dsub = dsub, pc = pc)
        }

      }

      fc = cStep(keep$dsub, K[i], keep$pc, keep$idxPC - 1L, keep$idxC - 1L)

      clustering = integer(N)
      clustering[keep$idx] = fc

      tempASW = .ASWCpp(clustering, dx, N, K[i])

      if(tempASW > bestASW){
        bestASW = tempASW
        bestClustering = clustering
      }

    }

    clusterings[, i] = bestClustering + 1L
    asw[i] = bestASW

  }

  .aswResult(clusterings, asw, K, "scalOSil", match.call(), list(n = n))

}
