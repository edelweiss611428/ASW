#' @name Init
#' @title Initialization methods for the Optimum Silhouette algorithm
#'
#' @description  This function computes an initialization for the Optimum Silhouette algorithm.
#'
#' @usage Init(dx, k, initMethod)
#'
#' @param dx  dx A "dist" object, which can be computed using stats::dist().
#' @param k An integer scalar specifying the number of clusters.
#' @param initMethod A character vector (or string) specifying initialization methods. Options include any
#' combination of "pam", "average", "single", "complete", "ward.D", "ward.D2", "mcquitty", "median", and "centroid".
#' By default, initMethod = "average".
#'
#' @return
#' \describe{
#' \item{clustering}{An initialized clustering.}
#' \item{asw}{The ASW associated with the initialized clustering.}
#' \item{method}{The "best" initialization method.}
#' }
#'
#' @details This function computes an initialization for the Optimum Silhouette algorithm, but it can be used as
#' a stand-alone clustering method.
#'
#'
#' @examples
#' x = faithful
#' dx = dist(x)
#' Initres = Init(dx, 2, c("pam", "average", "complete"))
#' plot(x, col = Initres$clustering, pch = 4)
#' print(paste(Initres$method, "achieves the highest ASW value"))
#'
#' @references
#' Batool, F. and Hennig, C., 2021. Clustering with the average silhouette width. Computational Statistics & Data Analysis, 158, p.107190.
#' Batool, F., 2019. Initialization methods for optimum average silhouette width clustering. arXiv preprint arXiv:1910.08644.
#'
#' @importFrom cluster pam
#' @importFrom stats hclust cutree dist
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

Init = function(dx, k, initMethod = "average"){

  if(inherits(dx, "dist") == TRUE){
    N = attr(dx, "Size")
  } else{
    stop("Init only inputs a distance matrix of class 'dist'")
  }

  if((!is.numeric(k)) | (length(k)!=1)){
    stop("k must be an integer")
  }

  k = as.integer(k)
  if(k > N){
    stop("k cannot be larger than the number of observations")
  } else if(k == 1){
    stop("k must be larger than 1")
  }

  if(!is.vector(initMethod, "character")){
    stop("initMethod must be a character vector (or string) specifying initialization methods!")
  }

  if(length(initMethod) == 0){
    stop("At least one initMethod needs to be specified!")
  }

  supportedMethods = c("pam", "average", "single", "complete", "ward.D",
                       "ward.D2", "mcquitty", "median", "centroid")

  if(!setequal(intersect(initMethod,supportedMethods), initMethod)){
    stop("initMethod contains unsupported methods!")
  }

  bestASW = -1

  for(i in seq_along(initMethod)){

    if(initMethod[i] == "pam"){
      tempClustering = pam(dx, k)$clustering
      tempASW = .ASWCpp(tempClustering-1L, dx, N, k)
    } else{
      tempClustering = cutree(hclust(dx, initMethod[i]),k)
      tempASW = .ASWCpp(tempClustering-1L, dx, N, k)
    }

    if(tempASW > bestASW){
      bestASW = tempASW
      bestClustering = tempClustering
      bestMethod = initMethod[i]
    }
  }

  return(list(clustering = bestClustering, asw = bestASW, method = bestMethod))

}
