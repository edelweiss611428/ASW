#' @name effOSil
#' @title The Efficient Optimum Silhouette algorithm
#'
#' @description  This function implements the Efficient Optimum Silhouette (effOSil) algorithm.
#'
#' @usage effOSil(dx, K, initMethod, variant)
#'
#' @param dx A "dist" object, which can be computed using stats::dist().
#' @param K An integer vector (or scalar) specifying the numbers of clusters. By default, K = 2:12.
#' @param initMethod A character vector (or string) specifying initialization methods.
#' By default, initMethod = "average". See ?Init for more details.
#' @param variant A character string specifying a variant. Options include "efficient" and "original".
#' If variant = "original", the original OSil algorithm is used. If variant = "efficient", effOSil is used.
#' By default, variant = "efficient".
#'
#' @return
#' \describe{
#' \item{best_clustering}{The effOSil clustering achieving the highest ASW value.}
#' \item{best_asw}{The highest ASW value.}
#' \item{k}{The estimated number of clusters.}
#' \item{clusterings}{The effOSil clustering solutions for all k in K.}
#' \item{asw}{The ASW values associated with the effOSil clusterings.}
#' \item{nIter}{The numbers of iterations needed for convergence.}
#' }
#'
#' @details
#' This function implements the Efficient Optimum Silhouette (effOSil) algorithm, an O(N) runtime improvement of
#' the original, computationally expensive OSil algorithm proposed by Batool & Hennig (2021) (N is the number of
#' observations). An implementation of the original OSil algorithm is also available for run time comparisions.
#'
#'
#' @examples
#' dx = dist(faithful)
#' effC = effOSil(dx, 2:8)
#' par(mfrow = c(2,1))
#' plot(faithful, col = effC$best_clustering, pch = 4)
#' plot(2:8, effC$asw, type = "l", xlab = "k", ylab = "ASW")
#'
#' @references
#' Batool, F. and Hennig, C., 2021. Clustering with the average silhouette width. Computational Statistics & Data Analysis, 158, p.107190.
#'
#' @importFrom stats dist
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

effOSil = function(dx, K = 2:12, initMethod = "average", variant = "efficient"){

  if(inherits(dx, "dist") == TRUE){
    N = attr(dx, "Size")
  } else{
    stop("effOSil only inputs a distance matrix of class 'dist'!")
  }

  nK = length(K)

  if((!is.numeric(K)) | (nK == 0)){
    stop("K must be an integer vector!")
  }

  K = as.integer(K)
  nuniqueK = length(unique(K))
  minK = min(K)
  maxK = max(K)

  if(nuniqueK != nK){
    stop("Duplicated number of clusters!")
  } else if(minK <= 1){
    stop("The number of clusters must be larger than 1!")
  } else if(maxK > N){
    stop("The number of clusters cannot be larger than the number of observations!")
  }

  if(length(variant) != 1){
    stop("Only ONE variant could be specified!")
  }

  clusterings = matrix(integer(N*nK), nrow = N)
  nIter = integer(nK)
  asw = numeric(nK)
  colnames(clusterings) = K
  names(asw) = K
  names(nIter) = K

  if(variant == "efficient"){
    FUNC = function(dx, initC, N, k){return(.effOSilCpp(dx, initC, N, k))}
  } else if(variant == "original"){
    FUNC = function(dx, initC, N, k){return(.OSilCpp(dx, initC, N, k))}
  } else{
    stop("The variant is not supported!")
  }


  for(i in 1:nK){
    init = Init(dx, K[i], initMethod)$clustering - 1L
    OSilres = FUNC(dx, init, N, K[i])
    clusterings[,i] = OSilres$Clustering
    asw[i] = OSilres$ASW
    nIter[i] = OSilres$nIter

  }

  idx_max = which.max(asw)
  best_asw = asw[idx_max]
  best_clustering = clusterings[,idx_max]
  k = K[idx_max]

  return(list(best_clustering = best_clustering, best_asw = best_asw, k = k,
              clusterings = clusterings, asw = asw, nIter = nIter))

}



