#' Clustered Dendrogram
#'
#' Generates dendrogram for the samples
#' @param nano The nano object
#' @param distMethod The method for distance calculation. Defaults to "euclidean".
#' @param clustMethod The clustering method. Defaults to "complete".
#' @keywords cluster dendrogram
#' @export
#' @examples
#' plotDendrogram()
plotDendrogram = function(nano, distMethod="euclidean", clustMethod="complete"){
  counts <- removeAnnot(nano$counts)
  distMat <- dist(t(counts), method = distMethod)
  hc <- hclust(distMat, method=clustMethod)
  dend <- as.dendrogram(hc)
  plot(dend, horiz=TRUE)
}

#
# # Experimental plotting using ggdendro and ggplot2
# dhc <- as.dendrogram(hc)
# ddata <- dendro_data(dhc, type="rectangle")
# p <- ggplot(segment(ddata)) +
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#   coord_flip() +
#   scale_y_reverse(expand = c(0.2, 0)) +
#   theme_dendro()
# p
