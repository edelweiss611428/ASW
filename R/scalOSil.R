#' @name scalOSil
#' @title The Scalable Optimum Silhouette algorithm
#'
#' @description  This function implements the Scalable Optimum Silhouette algorithm.
#'
#' @usage scalOSil(dx, K, n, ns, rep, initMethod, variant)
#'
#' @param dx A "dist" object, which can be computed using stats::dist().
#' @param K An integer vector (or scalar) specifying the numbers of clusters. By default, K = 2:12.
#' @param n An integer specifying the sample size. If not specified (NULL), n is set to 0.2*N where
#' N is the number of observations.
#' @param ns An integer specifying the number of random samples used in each instance. By default, ns = 1.
#' @param rep An integer specifying the number of scalOSil instances. By default, rep = 10.
#' @param initMethod A character vector (or string) specifying initialization methods. By default,
#' initMethod = "average". See ?Init for more details.
#' @param variant A character string specifying a variant. Options include "scalable" and "original".
#' If variant = "original", the original FOSil algorithm is used. If variant = "scalable", scalOSil is used.
#' By default, variant = "scalable".
#'
#' @return
#' \describe{
#' \item{best_clustering}{The scalOSil clustering achieving the highest ASW value.}
#' \item{best_asw}{The highest ASW value.}
#' \item{k}{The estimated number of clusters.}
#' \item{clusterings}{The scalOSil clustering solutions for all k in K.}
#' \item{asw}{The ASW values associated with the scalOSil clusterings.}
#' }
#'
#' @details
#' This function implements the Scalable Optimum Silhouette (scalOSil) algorithm, an O(n) runtime improvement of
#' the original, computationally expensive Fast OSil (FOSil) algorithm proposed by Batool & Hennig (2021) (n is
#' the sample size). An implementation of the original FOSil algorithm is also available for run time comparisions.
#'
#'
#' @examples
#' dx = dist(faithful)
#' scalC = scalOSil(dx, 2:8)
#' par(mfrow = c(2,1))
#' plot(faithful, col = scalC$best_clustering, pch = 4)
#' plot(2:8, scalC$asw, type = "l", xlab = "k", ylab = "ASW")
#'
#' @references
#' Batool, F. and Hennig, C., 2021. Clustering with the average silhouette width. Computational Statistics & Data Analysis, 158, p.107190.
#'
#' @importFrom stats dist
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

scalOSil = function(dx, K = 2:12, n = NULL, ns = 1, rep = 10, initMethod = "average", variant = "scalable"){

  if(inherits(dx, "dist") == TRUE){
    N = attr(dx, "Size")
  } else{
    stop("scalOSil only inputs a distance matrix of class 'dist'!")
  }

  if(is.null(n)){
    n = ceiling(0.2*N)
  } else{

    if((!is.numeric(n)) | (length(n)!=1)){
      stop("n must be an integer")
    }
    n = as.integer(n)
    if(n < 2){
      stop("n must be larger than 1")
    }
  }

  if((!is.numeric(rep)) | (length(rep)!=1)){
    stop("rep must be an integer")
  }
  rep = as.integer(rep)
  if(rep < 1){
    stop("rep must be larger than zero")
  }

  if((!is.numeric(ns)) | (length(ns)!=1)){
    stop("ns must be an integer")
  }
  ns = as.integer(ns)
  if(ns < 1){
    stop("ns must be larger than zero")
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
  } else if(maxK > n){
    stop("The number of clusters cannot be larger than the sample size!")
  }

  if(length(variant) != 1){
    stop("Only ONE variant should be specified!")
  }

  clusterings = matrix(integer(N*nK), nrow = N)
  asw = numeric(nK)
  colnames(clusterings) = K
  names(asw) = K

  if(variant == "scalable"){
    PC_Step = function(dx1, iC, n, k){return(.scalOSil_PC(dx1, iC, n, k))}
    C_Step = function(dx, dx1, k, PCres, idxPC, idxC, n, n2, N){
      return(.scalOSil_C(dx, k, PCres, idxPC, idxC, n, n2, N))
    }
  } else if(variant == "original"){
    PC_Step = function(dx1, iC, n, k){return(.FOSil_PC(dx1, iC, n, k))}
    C_Step = function(dx, dx1, k, PCres, idxPC, idxC, n, n2, N){
      return(.FOSil_C(dx, dx1, k, PCres, idxPC, idxC, n, n2, N))
    }
  } else{
    stop("The variant is not supported!")
  }

  for(i in 1:nK){


    bestASW = -1

    for(j in 1:rep){
      best_PCASW = -1
      for(l in 1:ns){
        temp_idxVec = sample.int(N)
        temp_idxPC = temp_idxVec[1:n]
        temp_idxC = temp_idxVec[(n+1):N]
        temp_dx1 = .subDistCpp(dx, temp_idxPC-1L, FALSE, FALSE, N, n)
        iC = Init(temp_dx1, K[i], initMethod)$clustering-1L
        temp_PC = PC_Step(temp_dx1, iC, n, K[i])

        if(temp_PC$ASW > best_PCASW){
          best_PCASW = temp_PC$ASW
          idxVec = temp_idxVec
          idxPC = temp_idxPC
          idxC = temp_idxC
          dx1 = temp_dx1
          PCres = temp_PC
        }

      }

      FCres = C_Step(dx, dx1, K[i], PCres, idxPC-1L, idxC-1L, n, N-n, N)
      FCres[idxVec] = FCres
      tempASW = .ASWCpp(FCres, dx, N,K[i])
      if(tempASW > bestASW){
        bestASW = tempASW
        bestClustering = FCres
      }
    }

    clusterings[,i] = FCres
    asw[i] = bestASW

  }

  idx_max = which.max(asw)
  best_asw = asw[idx_max]
  best_clustering = clusterings[,idx_max]
  k = K[idx_max]

  return(list(best_clustering = best_clustering+1L, best_asw = best_asw, k = k,
              clusterings = clusterings+1L, asw = asw))

}
