#########################
## Utility Functions ####
#########################
#' Geometric Mean
#'
#' Compute geometric mean
#' @param x Numeric data
#' @param na.rm Logic, if NAs should be removed or not. Defaults to TRUE
#' @keywords geometric mean
#' @export
#' @examples
#' geo.mean()
geo.mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#' Coefficient of variation
#'
#' Compute the coefficient of variation
#' @param vec Numeric vector
#' @keywords cv
#' @export
#' @examples
#' getCV()
getCV = function(vec){
  vec.sd <- sd(vec)
  vec.mean <- mean(vec)
  vec.cv <- vec.sd/vec.mean
  vec.cv
}

#' Remove annotation columns
#'
#' Removes the annotation columns from a nano count matrix
#' @param counts Counts from a nano object
#' @keywords utility
#' @export
#' @examples
#' removeAnnot()
removeAnnot = function(counts){
  rownames(counts) <- counts$Name
  counts[,4:ncol(counts)]
}

#' Remove counts below threshold
#'
#' Removes the annotation columns from a nano count matrix
#' @param counts Counts from a nano object
#' @param countCutoff The expression threshold
#' @keywords utility
#' @export
#' @examples
#' removeNonExpressed()
removeNonExpressed = function(counts, countCutoff=5){
  mean.counts <- apply(X = counts, MARGIN = 1, FUN = mean)
  ncounts <- counts[mean.counts >= countCutoff,]
}
