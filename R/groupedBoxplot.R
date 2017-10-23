#' Grouped Boxplot
#'
#' Generate Boxplot of all expressed counts, divided into groups
#' @param nano The nano object
#' @param groups A vector defining the groups of the samples. Default: all samples belong to the same group.
#' @param countCutoff The mean count value a gene must have to be considered for the plot.
#' @keywords boxplot
#' @export
#' @examples
#' plotGroupBoxplot()
plotGroupBoxplot = function(nano, groups = NA, title = "Title", countCutoff = 3) {
  x <- nano$counts[grepl("Endogenous", nano$counts$CodeClass),]
  x <- removeAnnot(x)
  # If no group information is given, generate random group
  if (is.na(groups)){
    groups <- rep("All", ncol(x))
  }
  # Remove genes without expression
  x.mean <- apply(X = x, MARGIN = 1, FUN = mean)
  x <- x[x.mean >= countCutoff,]
  # Pseudo counts
  x <- x + 1
  no.rows <- nrow(x)
  mdata <- melt(x)
  groups <- rep(groups, each=no.rows)
  mdata$group <- groups
  mdata <- mdata[with(mdata, order(mdata$group)),]

  p <- ggplot(mdata, aes(x=variable, y=value, color=group)) +
    geom_boxplot() +
    theme(axis.text.x=element_text(size=5, angle = 90, vjust = 0.5)) +
    theme(legend.position="top") +
    facet_grid(~group,scales="free", space="free") +
    scale_y_log10(name = "Log10 Counts") +
    labs(title = title) +
    xlab("Sample")
  p
}
