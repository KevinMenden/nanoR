###############
### Heatmap ###
###############
#' Expression Heatmap
#'
#' Plots clustered heatmap of expression values
#' Low expressed genes are removed
#' @param nano The nano object
#' @param countCutoff The mean expression cutoff to apply. Defaults to 3.
#' @keywords expression heatmap
#' @export
#' @examples
#' plotHeatmapExpression()
plotHeatmapExpression <- function(nano, groups, countCutoff = 3) {

  counts <- nano$counts
  counts <- counts[grepl("Endogenous", counts$CodeClass),]
  counts <- removeAnnot(counts)
  # Remove counts with too low expression
  mean.expression <- apply(X=counts, MARGIN = 1, FUN = mean)
  expressed <- mean.expression >= countCutoff
  counts <- counts[expressed,]
  log2.counts <- log2(counts + 1)

  annot <- data.frame(Group = groups)
  rownames(annot) <- colnames(log2.counts)

  pheatmap(log2.counts, scale="row", show_rownames = F, show_colnames = F, annotation_col = annot)

}
