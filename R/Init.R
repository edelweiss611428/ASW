#' @name Init
#' @title Initialisation for the Optimum Silhouette algorithms
#'
#' @description Runs one or more clustering methods and returns the solution
#' achieving the highest Average Silhouette Width (ASW), for use as an
#' initialisation by [effOSil] and [scalOSil].
#'
#' @param dx A `dist` object, as returned by [stats::dist()].
#' @param k An integer specifying the number of clusters.
#' @param initMethod A character vector of methods to try. Any combination of
#' `"pam"` and the agglomeration methods accepted by [stats::hclust()]:
#' `"average"`, `"single"`, `"complete"`, `"ward.D"`, `"ward.D2"`,
#' `"mcquitty"`, `"median"` and `"centroid"`. Defaults to `"average"`.
#'
#' @return A list with components:
#' \describe{
#'   \item{clustering}{An integer vector giving the clustering achieving the
#'     highest ASW among the methods tried.}
#'   \item{asw}{The ASW of that clustering.}
#'   \item{method}{The method that produced it.}
#'   \item{all_asw}{A named numeric vector of the ASW attained by each method.}
#' }
#'
#' @details Passing several methods gives a better starting point than any one
#' of them, at the cost of running each; Batool (2019) found this matters for
#' the quality of the final OSil solution. The function is also usable on its
#' own, as a way of picking among several clusterings by ASW.
#'
#' `"median"` and `"centroid"` assume squared Euclidean dissimilarities and can
#' produce inversions otherwise; [stats::hclust()] does not check this.
#'
#' @examples
#' dx = dist(scale(faithful))
#' fit = Init(dx, 2, c("pam", "average", "complete", "single"))
#'
#' fit$method
#' fit$all_asw
#' plot(faithful, col = fit$clustering, pch = fit$clustering)
#'
#' @references
#' Batool, F. (2019). Initialization methods for optimum average silhouette
#' width clustering. \emph{arXiv preprint} arXiv:1910.08644.
#' \doi{10.48550/arXiv.1910.08644}
#'
#' Batool, F. and Hennig, C. (2021). Clustering with the average silhouette
#' width. \emph{Computational Statistics & Data Analysis}, 158, 107190.
#' \doi{10.1016/j.csda.2021.107190}
#'
#' @seealso [effOSil], [scalOSil], [asw].
#'
#' @importFrom cluster pam
#' @importFrom stats hclust cutree setNames
#'
#' @author Minh Long Nguyen \email{edelweiss611428@@gmail.com}
#' @export

Init = function(dx, k, initMethod = "average"){

  N = .checkDist(dx)
  k = .checkCount(k, "k", lower = 2L)

  if(k > N){
    stop("The number of clusters cannot exceed the number of observations.", call. = FALSE)
  }

  initMethod = .checkInitMethod(initMethod)
  initMethod = unique(initMethod)

  supported = c("pam", "average", "single", "complete", "ward.D",
                "ward.D2", "mcquitty", "median", "centroid")
  unsupported = setdiff(initMethod, supported)

  if(length(unsupported) > 0L){
    stop(sprintf("Unsupported initialisation method(s): %s. Supported: %s.",
                 paste(unsupported, collapse = ", "),
                 paste(supported, collapse = ", ")), call. = FALSE)
  }

  allASW = setNames(numeric(length(initMethod)), initMethod)
  bestASW = -Inf
  bestClustering = NULL
  bestMethod = NA_character_

  for(i in seq_along(initMethod)){

    clustering = if(initMethod[i] == "pam"){
      pam(dx, k)$clustering
    } else {
      cutree(hclust(dx, method = initMethod[i]), k)
    }

    clustering = as.integer(clustering)
    tempASW = .ASWCpp(clustering - 1L, dx, N, k)
    allASW[i] = tempASW

    if(tempASW > bestASW){
      bestASW = tempASW
      bestClustering = clustering
      bestMethod = initMethod[i]
    }

  }

  list(clustering = bestClustering, asw = bestASW, method = bestMethod,
       all_asw = allASW)

}
