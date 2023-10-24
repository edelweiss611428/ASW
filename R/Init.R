#' @name Init
#' @title Initialization methods for Optimum Average Silhouette Width clustering algorithm
#'
#' @description  This function computes an initialized clustering for Optimum Average Silhouette Width clustering algorithms
#'
#' @usage Init(dx, k, initMethod = "average")
#'
#' @param dx  A "dist" object, which can be obtained by the "dist" function.
#' @param k The number of clusters.
#' @param initMethod A character vector specifying initialization methods. It must contain only supported methods:
#' one of the two combined methods "multiple1" and "multiple2"; or any combination of of "pam", "average", "single",
#' "complete", "ward.D", "ward.D2", "mcquitty", "median", and "centroid".
#'
#' @return
#' \describe{
#' \item{Clustering}{An initialized clustering.}
#' \item{ASW}{The ASW associated with the initialized clustering.}
#' \item{Method}{The "best" initialization method.}
#' }
#'
#' @details This function computes an initialized clustering for Optimum Average Silhouette Width clustering algorithms
#' by using the clustering methods specified by initMethod, and return the clustering maximize the ASW. The two combined methods
#' "multiple1" and "multiple2" are:
#' \describe{
#' \item{"multiple1"}{PAM and average linkage.}
#' \item{"multiple2"}{PAM, average linkage, single linkage, and the Ward's method (ward.D2).}
#' }
#'
#'
#' @examples
#' x = iris[,-5]
#' dx = dist(x)
#' Init(dx,3,"multiple1")
#'
#' @importFrom cluster pam
#' @importFrom stats hclust cutree dist
#'
#' @references
#' Batool, F. and Hennig, C., 2021. Clustering with the average silhouette width. Computational Statistics & Data Analysis, 158, p.107190.
#' Batool, F., 2019. Initialization methods for optimum average silhouette width clustering. arXiv preprint arXiv:1910.08644.
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
    stop("k must be a number")
  }

  k = as.integer(k)
  if(k > N){
    stop("k cannot be larger than the number of observations")
  } else if(k == 1){
    stop("k must be larger than 1")
  }

  if(!is.vector(initMethod, "character")){
    stop("initMethod must be a character vector of a string specifying initialzation methods!")
  }

  n_methods = length(initMethod)

  if(n_methods == 0){
    stop("At least one method needs to be specified!")
  }

  if(n_methods == 1){

    if(initMethod == "multiple1"){

      PAM_clustering = pam(dx, k)$clustering
      AL_clustering = cutree(hclust(dx,"average"),k)

      PAM_ASW = .ASWCpp(PAM_clustering-1L, dx, N, k)
      AL_ASW = .ASWCpp(AL_clustering-1L, dx, N, k)

      ind = which.max(c(PAM_ASW, AL_ASW))

      if(ind == 1){
        return(list(Clustering = PAM_clustering, ASW = PAM_ASW, Method = "pam"))
      } else{
        return(list(Clustering = AL_clustering, ASW = AL_ASW, Method = "Average-linkage"))
      }

    } else if (initMethod == "multiple2"){

      PAM_clustering = pam(dx, k)$clustering
      AL_clustering = cutree(hclust(dx,"average"),k)
      SL_clustering = cutree(hclust(dx,"single"),k)
      WM_clustering = cutree(hclust(dx,"ward.D2"),k)

      PAM_ASW = .ASWCpp(PAM_clustering-1L, dx, N, k)
      AL_ASW = .ASWCpp(AL_clustering-1L, dx, N, k)
      SL_ASW = .ASWCpp(SL_clustering-1L, dx, N, k)
      WM_ASW = .ASWCpp(WM_clustering-1L, dx, N, k)
      ind = which.max(c(PAM_ASW, AL_ASW, SL_ASW, WM_ASW))

      if(ind == 1){
        return(list(Clustering = PAM_clustering, ASW = PAM_ASW, Method = "pam"))
      } else if(ind == 2){
        return(list(Clustering = AL_clustering, ASW = AL_ASW, Method = "Average-linkage"))
      } else if(ind == 3){
        return(list(Clustering = SL_clustering, ASW = SL_ASW, Method = "Single-linkage"))
      } else{
        return(list(Clustering = WM_clustering, ASW = WM_ASW, Method = "Ward's method"))
      }

    }

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
  return(list(Clustering = bestClustering, ASW = bestASW, Method = bestMethod))

}
