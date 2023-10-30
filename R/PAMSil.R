#' @name PAMSil
#' @title The PAMSil algorithm
#'
#' @description  This function implements the PAMSil algorithm.
#'
#' @usage PAMSil(dx, K)
#'
#' @param dx A "dist" object, which can be computed using stats::dist().
#' @param K An integer vector (or scalar) specifying the numbers of clusters. By default, K = 2:12.
#'
#' @return
#' \describe{
#' \item{best_clustering}{The PAMSil clustering achieving the highest ASW value.}
#' \item{best_asw}{The highest ASW value.}
#' \item{best_medoids}{The medoids associated with the clustering maximize the ASW.}
#' \item{k}{The estimated number of clusters.}
#' \item{clusterings}{The PAMSil clustering solutions for all k in K.}
#' \item{asw}{The ASW values associated with the PAMSil clusterings.}
#' \item{medoids}{The medoids associated with the clustering solutions.}
#' \item{nIter}{The numbers of iterations needed for convergence.}
#' }
#'
#' @details
#' This function implements the PAMSil algorithm proposed by Van der Laan et al. (2003).
#' It is a k-medoids clustering algorithm whose objective function is the Average Silhouette Width.
#'
#'
#' @examples
#' dx = dist(faithful)
#' pamsilC = PAMSil(dx, 2:8)
#' par(mfrow = c(2,1))
#' plot(faithful, col = pamsilC$best_clustering, pch = 4)
#' plot(2:8, pamsilC$asw, type = "l", xlab = "k", ylab = "ASW")
#'
#' @references
#' Van der Laan, M., Pollard, K. and Bryan, J., 2003. A new partitioning around medoids algorithm. Journal of Statistical Computation and Simulation, 73(8), pp.575-584.
#'
#' @importFrom cluster pam
#' @importFrom stats dist
#'
#' @author Minh Long Nguyen \email{edelweiss611428@gmail.com}
#' @export

PAMSil = function(dx, K = 2:12){

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

  clusterings = matrix(integer(N*nK), nrow = N)
  nIter = integer(nK)
  asw = numeric(nK)
  medoids = vector("list", length = nK)
  colnames(clusterings) = K
  names(asw) = K
  names(nIter) = K
  names(medoids) = K

  for(i in 1:nK){

    PAM = pam(dx, K[i])
    PAMSilres = .PAMSilCpp(dx, PAM$clustering-1L, PAM$id.med-1L, N, K[i])
    clusterings[,i] = PAMSilres$Clustering
    asw[i] = PAMSilres$ASW
    nIter[i] = PAMSilres$nIter
    medoids[[i]] = PAMSilres$medoids

  }

  idx_max = which.max(asw)
  best_asw = asw[idx_max]
  best_clustering = clusterings[,idx_max]
  k = K[idx_max]
  best_medoids = medoids[[idx_max]]

  return(list(best_clustering = best_clustering, best_asw = best_asw, best_medoids = best_medoids, k = k,
              clusterings = clusterings, asw = asw, medoids = medoids, nIter = nIter))

}





