#' @name PAMSil
#' @title The PAMSil algorithm
#'
#' @description Clusters a distance matrix by maximising the Average Silhouette
#' Width (ASW) over partitions around medoids, using the PAMSil algorithm of
#' Van der Laan, Pollard & Bryan (2003).
#'
#' @param dx A `dist` object, as returned by [stats::dist()].
#' @param K An integer vector specifying the numbers of clusters to consider.
#' Defaults to `2:12`.
#'
#' @return An object of class `"ASW"`, a list with components:
#' \describe{
#'   \item{best_clustering}{An integer vector giving the clustering that
#'     achieves the highest ASW.}
#'   \item{best_asw}{The highest ASW value.}
#'   \item{best_medoids}{The medoids of the clustering maximising the ASW.}
#'   \item{k}{The estimated number of clusters.}
#'   \item{clusterings}{An integer matrix with one column per element of `K`.}
#'   \item{asw}{A numeric vector of ASW values, one per element of `K`.}
#'   \item{medoids}{A list of medoid index vectors, one per element of `K`.}
#'   \item{nIter}{An integer vector of iteration counts to convergence.}
#'   \item{method}{The name of the algorithm used.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details PAMSil is a k-medoid algorithm whose objective function is the ASW.
#' It is initialised with [cluster::pam()] and then repeatedly replaces the
#' medoid whose replacement increases the ASW most, stopping when no replacement
#' improves it. Because it searches over medoid assignments rather than over
#' arbitrary partitions, it is more constrained than [effOSil] and will usually
#' attain a lower ASW.
#'
#' @examples
#' x = scale(faithful)
#' dx = dist(x)
#' fit = PAMSil(dx = dx, K = 2:6)
#'
#' oldpar = par(mfrow = c(1, 2))
#' plot(faithful, col = fit$best_clustering, pch = fit$best_clustering)
#' plot(2:8, fit$asw, type = "l", xlab = "k", ylab = "ASW")
#' par(oldpar)
#'
#' @references
#' Van der Laan, M., Pollard, K. and Bryan, J. (2003). A new partitioning around
#' medoids algorithm. \emph{Journal of Statistical Computation and Simulation},
#' 73(8), 575-584. \doi{10.1080/0094965031000136012}
#'
#' @seealso [effOSil], [scalOSil], [asw], [Silhouette].
#'
#' @author Minh Long Nguyen \email{edelweiss611428@@gmail.com}
#' @export

PAMSil = function(dx, K = 2:12){

  N = .checkDist(dx)
  K = .checkK(K, N)

  nK = length(K)
  clusterings = matrix(NA_integer_, nrow = N, ncol = nK, dimnames = list(NULL, K))
  asw = setNames(numeric(nK), K)
  nIter = setNames(integer(nK), K)
  medoids = setNames(vector("list", nK), K)

  for(i in seq_len(nK)){
    init = pam(dx, K[i])
    res = .PAMSilCpp(dx, init$clustering - 1L, init$id.med - 1L, N, K[i])
    clusterings[, i] = res$Clustering
    asw[i] = res$ASW
    nIter[i] = res$nIter
    medoids[[i]] = res$medoids
  }

  idxMax = which.max(asw)

  .aswResult(clusterings, asw, K, "PAMSil", match.call(),
             list(best_medoids = medoids[[idxMax]], medoids = medoids, nIter = nIter))

}
