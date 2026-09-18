#' @name effOSil
#' @title The Efficient Optimum Silhouette algorithm
#'
#' @description Clusters a distance matrix by maximising the Average Silhouette
#' Width (ASW), using the Efficient Optimum Silhouette (effOSil) algorithm.
#'
#' @param dx A `dist` object, as returned by [stats::dist()].
#' @param K An integer vector specifying the numbers of clusters to consider.
#' Defaults to `2:12`.
#' @param initMethod A character vector of initialisation methods. Defaults to
#' `"average"`; to obtain a better initialisation in terms of the ASW, several
#' methods may be given, for instance
#' `c("single", "average", "complete", "pam")`. See [Init] for the available
#' methods.
#' @param variant Either `"efficient"` (the default), which uses effOSil, or
#' `"original"`, which uses the original OSil algorithm. Both return the same
#' clustering; `"original"` is provided for comparison and is far slower.
#'
#' @return An object of class `"ASW"`, a list with components:
#' \describe{
#'   \item{best_clustering}{An integer vector giving the clustering that
#'     achieves the highest ASW.}
#'   \item{best_asw}{The highest ASW value.}
#'   \item{k}{The estimated number of clusters.}
#'   \item{clusterings}{An integer matrix with one column per element of `K`.}
#'   \item{asw}{A numeric vector of ASW values, one per element of `K`.}
#'   \item{nIter}{An integer vector of iteration counts to convergence.}
#'   \item{method}{The name of the algorithm used.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details effOSil returns exactly the same clustering as the OSil algorithm of
#' Batool & Hennig (2021), but evaluates each candidate swap in \eqn{O(N)}
#' rather than \eqn{O(N^2)} time by caching, for every observation, the mean
#' distance to its own cluster and to the three nearest other clusters. This
#' gives an \eqn{O(N)} reduction in overall runtime, where \eqn{N} is the number
#' of observations.
#'
#' Both variants use steepest ascent: each iteration evaluates every
#' single-observation reassignment and applies the one that increases the ASW
#' most, stopping when no reassignment improves it. Clusters are never allowed
#' to become empty.
#'
#' On data containing duplicated observations the ASW can have exact ties, in
#' which case `"efficient"` and `"original"` may resolve a tied reassignment
#' differently and converge to different local optima of equal or near-equal
#' ASW.
#'
#' @examples
#' x = scale(faithful)
#' dx = dist(x)
#' fit = effOSil(dx = dx, K = 2:8)
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
#' @seealso [scalOSil] for larger data sets, [PAMSil], [Init], [asw], [Silhouette].
#'
#' @author Minh Long Nguyen \email{edelweiss611428@@gmail.com}
#' @export

effOSil = function(dx, K = 2:12, initMethod = "average",
                   variant = c("efficient", "original")){

  N = .checkDist(dx)
  K = .checkK(K, N)
  initMethod = .checkInitMethod(initMethod)
  variant = match.arg(variant)

  nK = length(K)
  clusterings = matrix(NA_integer_, nrow = N, ncol = nK, dimnames = list(NULL, K))
  asw = setNames(numeric(nK), K)
  nIter = setNames(integer(nK), K)

  osil = if(variant == "efficient") .effOSilCpp else .OSilCpp

  for(i in seq_len(nK)){
    init = Init(dx, K[i], initMethod)$clustering - 1L
    res = osil(dx, init, N, K[i])
    clusterings[, i] = res$Clustering
    asw[i] = res$ASW
    nIter[i] = res$nIter
  }

  .aswResult(clusterings, asw, K, "effOSil", match.call(), list(nIter = nIter))

}
