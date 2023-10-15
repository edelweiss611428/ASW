#' @name ASW
#' @title The Average Silhouette Width
#'
#' @description  This function computes the Average Silhouette Width (ASW).
#'
#' @usage ASW(C, dist)
#'
#' @param C A clustering solution. It must be an integer vector of k unique values 1,2,...,k.
#' @param dist  A "dist" object, which can be obtained by the "dist" function.
#'
#' @return The ASW of the clustering C with respect to the distance matrix dist.
#'
#' @references
#' Rousseeuw, P.J. (1987) Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. J. Comput. Appl. Math., 20, 53–65.
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

ASW = function(C, dist){

  if(inherits(dist, "dist") == TRUE){
    N = (1+sqrt(1+8*length(dist)))/2
  } else{
    stop("ASW only inputs a distance matrix of class 'dist'.")
  }

  k = length(unique(C))
  C = as.integer(C)

  if(length(C) != N | min(C) != 1 | max(C) != k){
    stop("Not a valid clustering!")
  }

  return(.ASWCpp(C-1L, dist, N, k))

}


