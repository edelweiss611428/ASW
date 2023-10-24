#' @name scalOSil
#' @title The Scalable Optimum Silhouette algorithm
#'
#' @description  This function implements the Scalable Optimum Silhouette (scalOSil) algorithm.
#'
#' @usage scalOSil(dx, k, n = "default", rep = 5, initMethod = "average")
#'
#' @param dx  A "dist" object, which can be obtained by the "dist" function.
#' @param k The number of clusters.
#' @param n The sample size. By default, n = ceiling(0.2*N).
#' @param rep The number of scalOSil instances. By default, rep = 5.
#' @param initMethod A character vector specifying initialization methods. It must contain only supported methods:
#' one of the two combined methods "multiple1" and "multiple2"; or any combination of of "pam", "average", "single",
#' "complete", "ward.D", "ward.D2", "mcquitty", "median", and "centroid". See ?Init for more details.
#'
#' @return
#' \describe{
#' \item{Clustering}{Final clustering.}
#' \item{ASW}{The ASW of the scalOSil clustering w.r.t. dx.}
#' }
#'
#' @details
#' The scalOSil algorithm is an approximation algorithm of effOSil based on subsetting. It is an improved version of FOSil (Batool & Hennig 2021).
#' Both the algorithms consists of two steps: partial clustering (PC) and classification (C).
#'
#' In the PC-step of scalOSil, effOSil is applied to a random subset of the dataset, obtaining a subset S and its effOSil
#' clustering C_S. In the C-step of scalOSil, each unassigned data point is classified into one of the clusters in C_S in such
#' a way that the ASW is maximized.
#'
#' Unlike FOSil, scalOSil runs many instances, specified by the parameter "rep", and for each instance, scalOSil only runs the
#' PC-step once. Moreover, the PC-step of scalOSil scales quadratically in n, while that of FOSil scales cubically in n, and
#' the C-step of scalOSil scales linearly in n, while that of FOSil scales quadratically in n. These allow scalOSil to handle
#' much larger datasets.
#'
#' @examples
#' x = iris[,-5]
#' dx = dist(x)
#' scalOSil_clustering = scalOSil(dx, 3, initMethod = "average")
#' plot(x, col = scalOSil_clustering$Clustering)
#'
#' @references
#' Batool, F. and Hennig, C., 2021. Clustering with the average silhouette width. Computational Statistics & Data Analysis, 158, p.107190.
#'
#' @importFrom cluster pam
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

scalOSil = function(dx, k, n = "default", rep = 5, initMethod = "average"){

  if(inherits(dx, "dist") == TRUE){
    N = attr(dx, "Size")
  } else{
    stop("scalOSil only inputs a distance matrix of class 'dist'")
  }

  if(n == "default"){
    n = ceiling(0.2*N)
  } else{
    if((!is.numeric(n)) | (length(n)!= 1)){
      stop("If n is not 'default', it must be a number!")
    }
  }

  if((!is.numeric(k)) | (length(k)!=1)){
    stop("k must be a number")
  }

  k = as.integer(k)
  if(k > n){
    stop("k cannot be larger than n")
  } else if(k == 1){
    stop("k must be larger than 1")
  }

  if((!is.numeric(rep)) | (length(rep)!=1)){
    stop("rep must be a number")
  }
  rep = as.integer(rep)

  bestASW = -1

  for(i in 1:rep){
    idxVec = sample.int(N)
    idxPC = idxVec[1:n]
    idxC = idxVec[(n+1):N]
    dx1 = .subDistCpp(dx, idxPC, FALSE, FALSE, N = N, n = n)
    iC = Init(dx1, k, initMethod)$Clustering -1L
    PCres = .scalOSil_PC(dx1, iC, n, k)
    FCres = .scalOSil_C(dx,k, PCres, idxPC-1L, idxC-1L, n, N-n, N)
    FCres[idxVec] = FCres
    tempASW = .ASWCpp(FCres, dx, N,k)

    if(tempASW > bestASW){
      bestASW = tempASW
      bestClustering = FCres
    }
  }

  return(list(Clustering = bestClustering+1, ASW = bestASW))
}

