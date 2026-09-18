#' @name asw
#' @title Average Silhouette Width of a clustering
#'
#' @description Computes the Average Silhouette Width (ASW) of a clustering.
#'
#' @param clustering A vector of cluster labels, one per observation. Any
#' labels may be used; they are matched to contiguous integers internally.
#' @param dx A `dist` object, as returned by [stats::dist()].
#'
#' @return A single numeric value, the ASW.
#'
#' @details Observations in singleton clusters are assigned a silhouette width
#' of 0, following Rousseeuw (1987). Coincident observations, for which both the
#' within- and between-cluster mean distances are zero, are also assigned 0.
#'
#' @examples
#' dx = dist(scale(faithful))
#' fit = effOSil(dx, K = 2:4)
#' asw(fit$best_clustering, dx)
#'
#' @references
#' Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the interpretation
#' and validation of cluster analysis. \emph{Journal of Computational and
#' Applied Mathematics}, 20, 53-65. \doi{10.1016/0377-0427(87)90125-7}
#'
#' @seealso [Silhouette] for the per-observation widths, [effOSil].
#'
#' @author Minh Long Nguyen \email{edelweiss611428@@gmail.com}
#' @export

asw = function(clustering, dx){

  N = .checkDist(dx)
  z = .zeroBased(clustering, N)
  .ASWCpp(z$labels, dx, N, z$k)

}

#' @name print.ASW
#' @title Print a clustering produced by the ASW package
#'
#' @param x An object of class `"ASW"`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly. Called for its side effect of printing.
#'
#' @examples
#' effOSil(dist(scale(faithful)), K = 2:4)
#'
#' @export

print.ASW = function(x, ...){

  cat(sprintf("%s clustering\n\n", x$method))
  cat(sprintf("  observations : %d\n", nrow(x$clusterings)))
  cat(sprintf("  k considered : %s\n", paste(colnames(x$clusterings), collapse = ", ")))
  cat(sprintf("  k selected   : %d\n", x$k))
  cat(sprintf("  ASW          : %.4f\n\n", x$best_asw))
  cat("  cluster sizes:\n")
  print(table(x$best_clustering))

  invisible(x)

}
