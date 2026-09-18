#' Validate a distance matrix
#'
#' @param dx An object to validate.
#' @return The number of observations, as an integer.
#' @noRd

.checkDist = function(dx){

  if(!inherits(dx, "dist")){
    stop("`dx` must be a 'dist' object, as returned by stats::dist().", call. = FALSE)
  }
  if(anyNA(dx)){
    stop("`dx` contains missing values.", call. = FALSE)
  }
  if(min(dx) < 0){
    stop("`dx` contains negative dissimilarities.", call. = FALSE)
  }

  as.integer(attr(dx, "Size"))

}

#' Validate a vector of cluster numbers
#'
#' @param K An object to validate.
#' @param upper The largest admissible number of clusters.
#' @param what A label for `upper`, used in the error message.
#' @return `K`, as an integer vector.
#' @noRd

.checkK = function(K, upper, what = "the number of observations"){

  if(!is.numeric(K) || length(K) == 0L || anyNA(K)){
    stop("`K` must be a non-empty numeric vector without missing values.", call. = FALSE)
  }

  K = as.integer(K)

  if(anyDuplicated(K)){
    stop("`K` must not contain duplicated values.", call. = FALSE)
  }
  if(min(K) < 2L){
    stop("The number of clusters must be at least 2.", call. = FALSE)
  }
  if(max(K) > upper){
    stop(sprintf("The number of clusters cannot exceed %s (%d).", what, upper), call. = FALSE)
  }

  K

}

#' Validate a single positive integer
#'
#' @param x An object to validate.
#' @param name The argument name, used in the error message.
#' @param lower The smallest admissible value.
#' @return `x`, as a single integer.
#' @noRd

.checkCount = function(x, name, lower = 1L){

  if(!is.numeric(x) || length(x) != 1L || is.na(x)){
    stop(sprintf("`%s` must be a single integer.", name), call. = FALSE)
  }

  x = as.integer(x)

  if(x < lower){
    stop(sprintf("`%s` must be at least %d.", name, lower), call. = FALSE)
  }

  x

}

#' Validate initialisation methods
#'
#' @param initMethod An object to validate.
#' @return `initMethod`, unchanged.
#' @noRd

.checkInitMethod = function(initMethod){

  if(!is.character(initMethod) || length(initMethod) == 0L || anyNA(initMethod)){
    stop("`initMethod` must be a non-empty character vector.", call. = FALSE)
  }

  initMethod

}

#' Map an arbitrary clustering to zero-based contiguous labels
#'
#' @param clustering A vector of cluster labels.
#' @param N The expected number of observations.
#' @return A list with the zero-based labels and the number of clusters.
#' @noRd

.zeroBased = function(clustering, N){

  if(length(clustering) != N){
    stop(sprintf("`clustering` must have length %d.", N), call. = FALSE)
  }
  if(anyNA(clustering)){
    stop("`clustering` contains missing values.", call. = FALSE)
  }

  lvls = sort(unique(clustering))

  if(length(lvls) < 2L){
    stop("`clustering` must contain at least 2 clusters.", call. = FALSE)
  }

  list(labels = match(clustering, lvls) - 1L, k = length(lvls))

}

#' Assemble the return value of a clustering routine
#'
#' @param clusterings A matrix of clusterings, one column per k.
#' @param asw The ASW values associated with the columns.
#' @param K The numbers of clusters.
#' @param method The name of the algorithm.
#' @param call The matched call.
#' @param extra Further named elements to append.
#' @return A list of class `"ASW"`.
#' @noRd

.aswResult = function(clusterings, asw, K, method, call, extra = list()){

  idxMax = which.max(asw)

  out = c(list(best_clustering = clusterings[, idxMax],
               best_asw = unname(asw[idxMax]),
               k = K[idxMax],
               clusterings = clusterings,
               asw = asw),
          extra,
          list(method = method, call = call))

  structure(out, class = "ASW")

}
