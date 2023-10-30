#' @name Silhouette
#' @title Silhouette Width
#'
#' @description  This function computes the Silhouette Widths for all data points in the dataset.
#'
#' @usage Silhouette(C, dx)
#'
#' @param C An integer vector specifying a k-partition of the dataset. min(C) must be 1
#' and max(C) must be k.
#' @param dx  A "dist" object, which can be computed using stats::dist().
#'
#' @return A numeric matrix of class "silhouette" containing three columns
#' \describe{
#' \item{cluster}{A clustering of the dataset.}
#' \item{neighbor}{The clustering labels of the nearest clusters for all data points.}
#' \item{sil_width}{The silhouette widths of data points.}
#' }
#'
#'
#' @examples
#' library("cluster")
#' dx = dist(faithful)
#' C = pam(dx, 2)$clustering
#' plot(Silhouette(C,dx))
#'
#' @references
#' Rousseeuw, P.J. (1987) Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. J. Comput. Appl. Math., 20, 53–65.
#'
#' @importFrom cluster pam
#' @importFrom stats dist
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

Silhouette = function (C, dx){
  cll = match.call()
  if (inherits(dx, "dist") == TRUE) {
    N = attr(dx, "Size")
  }
  else {
    stop("Silhouette only inputs a distance matrix of class 'dist'.")
  }
  k = length(unique(C))
  C = as.integer(C)
  if (length(C) != N | min(C) != 1 | max(C) != k) {
    stop("Not a valid clustering!")
  }
  SW = .SWCpp(C - 1L, dx, N, k)

  wds = cbind(cluster = C, neighbor = (SW$neighbor+1), sil_width = SW$sil_width)
  attr(wds, "Ordered") <- FALSE
  attr(wds, "call") <- cll
  class(wds) <- "silhouette"
  return(wds)
}

