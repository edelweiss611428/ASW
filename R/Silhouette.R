#' @name Silhouette
#' @title Silhouette widths
#'
#' @description Computes the silhouette width of every observation, together
#' with its neighbouring cluster.
#'
#' @param C An integer vector giving a clustering, with labels `1` to `k`.
#' @param dx A `dist` object, as returned by [stats::dist()].
#'
#' @return A numeric matrix of class `"silhouette"` with one row per
#' observation and three columns:
#' \describe{
#'   \item{cluster}{The observation's cluster.}
#'   \item{neighbor}{The nearest other cluster.}
#'   \item{sil_width}{The silhouette width.}
#' }
#' The result can be passed to the `plot()`, `print()` and `summary()` methods
#' of the \pkg{cluster} package.
#'
#' @details Observations in singleton clusters are given a silhouette width of
#' 0, following Rousseeuw (1987); their `neighbor` is still reported, since it
#' remains well defined. Coincident observations, for which both the within- and
#' between-cluster mean distances are zero, are also given a width of 0.
#'
#' @examples
#' dx = dist(scale(faithful))
#' fit = effOSil(dx, K = 2:5)
#' sw = Silhouette(fit$best_clustering, dx)
#'
#' summary(sw)
#' plot(sw)
#'
#' @references
#' Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the interpretation
#' and validation of cluster analysis. \emph{Journal of Computational and
#' Applied Mathematics}, 20, 53-65. \doi{10.1016/0377-0427(87)90125-7}
#'
#' @seealso [asw] for the average, [effOSil], [cluster::silhouette()].
#'
#' @author Minh Long Nguyen \email{edelweiss611428@@gmail.com}
#' @export

Silhouette = function(C, dx){

  cll = match.call()
  N = .checkDist(dx)

  if(length(C) != N){
    stop(sprintf("`C` must have length %d.", N), call. = FALSE)
  }
  if(!is.numeric(C) || anyNA(C)){
    stop("`C` must be an integer vector without missing values.", call. = FALSE)
  }

  C = as.integer(C)
  k = length(unique(C))

  if(min(C) != 1L || max(C) != k){
    stop("`C` must use the labels 1 to k, each at least once.", call. = FALSE)
  }

  SW = .SWCpp(C - 1L, dx, N, k)

  wds = cbind(cluster = C, neighbor = SW$neighbor, sil_width = SW$sil_width)
  attr(wds, "Ordered") = FALSE
  attr(wds, "call") = cll
  class(wds) = "silhouette"

  wds

}
