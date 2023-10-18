#' @name effOSil
#' @title Efficient Optimum Silhouette Clustering Algorithm
#'
#' @description  effOSil implements the exact Optimum Silhouette (OSil) algorithm. It is
#' O(N) times faster than the original OSil algorithm at the cost of storing O(N) additional
#' values.
#'
#' @usage effOSil(dist, k, initClustering = NULL, initMethod = "average")
#'
#' @param dist  A "dist" object, which can be obtained by the "dist" function.
#' @param k The number of clusters.
#' @param initClustering A user-specified initialized clustering. It must be an integer vector of
#' k unique values 1,2,...,k. By defaualt, initClustering is set to NULL. If initClustering is NULL,
#' initMethod is used; otherwise, initClustering is used.
#' @param initMethod An initialization method. By default, effOSil uses average-linkage clustering
#' ("average") for initialization. Other options include well-established clustering algorithms
#' such as Partition Around Medoids ("pam"), single-linkage clustering ("single"), complete-linkage
#' clustering ("complete"), the Ward's method ("ward.D"), etc.
#'
#' @return
#' \describe{
#' \item{Clustering}{The OSil clustering solution.}
#' \item{ASW}{The ASW associated with the OSil clustering.}
#' \item{nIter}{The number of iterations needed for convergence.}
#' }
#'
#' @importFrom cluster pam
#' @importFrom stats hclust cutree
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

effOSil = function(dist, k, initClustering = NULL, initMethod = "average"){

  if(inherits(dist, "dist") == TRUE){
    N = attr(dist, "Size")
  } else{
    stop("effOSil only inputs a distance matrix of class 'dist'")
  }

  if(k > N){
    stop("The number of clusters cannot be larger than the
         number of observations")
  } else if(k == 1){
    stop("The number of clusters must be larger than 1")
  }

  if(!is.null(initClustering)){

    initClustering = as.integer(initClustering)

    if(length(initClustering) != N|
       length(unique(initClustering)) != k |
       min(initClustering) != 1 |
       max(initClustering) != k){
      stop("Not a valid initialized clustering!")
    }

    return(.effOSilCpp(dist, initClustering-1L, N,k))

  } else{

    if(!is.element(initMethod, c("pam", "average", "single", "complete",
                                 "ward.D", "ward.D2", "mcquitty", "median", "centroid"))){
      stop("The initMethod is not supported!")
    }

    if(initMethod == "pam"){
      initClustering = pam(dist, k)$clustering-1L
    } else{
      initClustering = cutree(hclust(dist, initMethod), k)-1L
    }

    return(.effOSilCpp(dist, initClustering, N,k))

  }
}


