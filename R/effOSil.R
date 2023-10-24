#' @name effOSil
#' @title The Efficient Optimum Silhouette algorithm
#'
#' @description  This function implements the exact Optimum Silhouette (OSil) algorithm.
#'
#' @usage effOSil(dx, k, initClustering = NULL, initMethod = "average")
#'
#' @param dx  A "dist" object, which can be obtained by the "dist" function.
#' @param k The number of clusters.
#' @param initClustering An initialized clustering. It must be an numeric vector of k unique values 1,2,...,k.
#' By default, initClustering is set to NULL. If initClustering is NULL, initMethod is used instead; otherwise, initClustering is used.
#' @param initMethod A character vector specifying initialization methods. It must contain only supported methods:
#' one of the two combined methods "multiple1" and "multiple2"; or any combination of of "pam", "average", "single",
#' "complete", "ward.D", "ward.D2", "mcquitty", "median", and "centroid". See ?Init for more details.
#'
#' @return
#' \describe{
#' \item{Clustering}{The OSil clustering solution.}
#' \item{ASW}{The ASW associated with the OSil clustering.}
#' \item{nIter}{The number of iterations needed for convergence.}
#' }
#'
#' @details
#' This function implements the exact Optimum Silhouette (OSil) algorithm proposed by Batool & Hennig (2021).
#' However, it is O(N) times faster than the original OSil algorithm at the cost of storing O(kN) additional values.
#'
#'
#' @examples
#' x = iris[,-5]
#' dx = dist(x)
#' effOSil_clustering = effOSil(dx, 3, initMethod = "average")
#' plot(x, col = effOSil_clustering$Clustering)
#'
#' @references
#' Batool, F. and Hennig, C., 2021. Clustering with the average silhouette width. Computational Statistics & Data Analysis, 158, p.107190.
#'
#' @importFrom cluster pam
#' @importFrom stats hclust cutree dist
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

effOSil = function(dx, k, initClustering = NULL, initMethod = "average"){

  if(inherits(dx, "dist") == TRUE){
    N = attr(dx, "Size")
  } else{
    stop("effOSil only inputs a distance matrix of class 'dist'")
  }

  if((!is.numeric(k)) | (length(k)!=1)){
    stop("k must be a number")
  }

  k = as.integer(k)
  if(k > N){
    stop("k cannot be larger than the number of observations")
  } else if(k == 1){
    stop("k must be larger than 1")
  }

  if(!is.null(initClustering)){

    initClustering = as.integer(initClustering)

    if(length(initClustering) != N|
       length(unique(initClustering)) != k |
       min(initClustering) != 1 |
       max(initClustering) != k){
      stop("Not a valid initialized clustering!")
    }

    return(.effOSilCpp(dx, initClustering-1L, N,k))

  } else{

    initClustering = Init(dx,k, initMethod)$Clustering - 1L
    return(.effOSilCpp(dx, initClustering, N,k))

  }
}


