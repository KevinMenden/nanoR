#' Volcano Plot
#'
#' Generates Volcano Plots for every contrast given
#' @param nano The nano object
#' @param groups Group information for the columns of the count matrix
#' @keywords principal component analysis pca
#' @export
#' @examples
#' plotPCA()
volcanoPlot = function(dge.result, contrast){
  dge.result$threshold <- as.factor(abs(dge.result$logFC) > 1 & dge.result$adj.P.Val < 0.05)
  dge.result$comparison <- contrast
  dge.result$geneID <- rownames(dge.result)

  volcano.plot <- ggplot(data=dge.result, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
    scale_x_continuous(limits = c(-12, 12)) + scale_y_continuous(limits = c(0, 20)) +
    geom_point(alpha=0.4, size=3) + xlab("log2 fold change") + ylab("-log10 P.Value") + ggtitle(contrast) +
    scale_colour_manual(values=c("#00B0E5","#FF7373"))

  return(volcano.plot)
  #ggsave(volcano.plot,filename=paste("Volcano_",cont,".pdf",sep=""))
}
