#' @name FOSil
#' @title The Fast Optimum Silhouette algorithm
#'
#' @description  This function implements the Fast Optimum Silhouette (FOSil) algorithm.
#'
#' @usage FOSil(dx, k, n = "default", ns = 25, initMethod = "average")
#'
#' @param dx  A "dist" object, which can be obtained by the "dist" function.
#' @param k The number of clusters.
#' @param n The sample size. By default, n = ceiling(0.2*N).
#' @param ns The number of scalOSil instances. By default, ns = 25.
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
#' The Fast Optimum Silhouette algorithm (FOSil; Batool & Hennig (2021)) is an approximation algorithm of OSil,
#' based on subsetting. It consists of two steps: partial clustering (PC) and classification (C).
#'
#' In the PC-step of FOSil, FOSil is applied to various subsets of equal size, in which the subset S and its OSil
#' clustering C_S maximizing the ASW is selected. In the C-step of FOSil, each unassigned data point is classified
#' into one of the clusters in C_S in such a way that the ASW is maximized.
#'
#' However, FOSil is still a computationally expensive algorithm. The PC-step of FOSil scales cubically in
#' n and the C-step of FOSil scales quadratically in n. scalOSil is an improved version of FOSil, which improves both
#' steps of FOSil by O(n) time, allowing us to handle much larger datasets (see ?scalOSil for more details).
#'
#' @examples
#' x = iris[,-5]
#' dx = dist(x)
#' FOSil_clustering = FOSil(dx, 3, initMethod = "average")
#' plot(x, col = FOSil_clustering$Clustering)
#'
#' @references
#' Batool, F. and Hennig, C., 2021. Clustering with the average silhouette width. Computational Statistics & Data Analysis, 158, p.107190.
#'
#' @importFrom cluster pam
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

FOSil = function(dx, k, n = "default", ns = 25, initMethod = "average"){

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

  if((!is.numeric(ns)) | (length(ns)!=1)){
    stop("rep must be a number")
  }
  ns = as.integer(ns)

  bestPCASW = -1

  for(i in 1:ns){
    tempidxVec = sample.int(N)
    tempidxPC = tempidxVec[1:n]
    tempidxC = tempidxVec[(n+1):N]
    tempdx1 = .subDistCpp(dx, tempidxPC-1L, FALSE, FALSE, N, n)

    iC = Init(tempdx1, k, initMethod)$Clustering-1

    tempPC = .FOSil_PC(tempdx1, iC, n, k)

    if(tempPC$ASW > bestPCASW){
      bestPCASW = tempPC$ASW
      idxVec = tempidxVec
      idxPC = tempidxPC
      idxC = tempidxC
      dx1 = tempdx1
      PCres = tempPC
    }

  }

  FCres = .FOSil_C(dx, dx1, k, PCres, idxPC-1L, idxC-1L, n, N-n, N)
  FCres[idxVec] = FCres
  ASWFC = .ASWCpp(FCres, dx, N, k)

  return(list(Clustering = FCres+1, ASW = ASWFC))
}

