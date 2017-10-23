########################
### PCA related code ###
########################
# pca plot, pca heatmap

# Find Hulls
find_hulls <- function(pca.df){
  df <- data.frame()
  gs <- levels(factor(pca.df$group))
  for (g in gs){
    g.df <- pca.df[pca.df$group == g,]
    tmp <- g.df[chull(g.df$PC1, g.df$PC2),]
    df <- rbind(df, tmp)
  }
  return(df)
}

# Create PCA dataframe
createPCAdf <- function(counts, groups){
  # Calculate principal components
  prc <- prcomp(t(counts),scale=T)
  pca.matrix <- prc$x
  pca.matrix <- cbind(as.character(groups),pca.matrix)
  colnames(pca.matrix)<- c("group",paste("PC",seq(1:dim(prc$rotation)[2]),sep=""))

  # Make PCA dataframe
  pca.df <- as.data.frame(pca.matrix)
  pca.df[,2:dim(pca.df)[2]] <- sapply(pca.df[,2:dim(pca.df)[2]], as.character)
  pca.df[,2:dim(pca.df)[2]] <- sapply(pca.df[,2:dim(pca.df)[2]], as.numeric)

  return(pca.df)
}


#' Principal Component Analysis Plots
#'
#' Generates PCA plots
#' @param nano The nano object
#' @param groups Group information for the columns of the count matrix
#' @param countCutoff Mean cut value used to cut off unexpressed genes. Defaults to 5
#' @keywords principal component analysis pca
#' @export
#' @examples
#' plotPCA()
plotPCA = function(nano, groups, countCutoff = 5){
  # Extract expressed counts
  counts <- nano$counts
  no.cols <- ncol(counts)
  counts <- counts[,4:no.cols]
  mean.counts <- apply(X = counts, MARGIN = 1, FUN = mean)
  counts <- counts[mean.counts >= countCutoff,]
  Group <- groups

  # Calculate principal components
  pca.df <- createPCAdf(counts, groups)


  # Plot pcas
  hulls <- find_hulls(pca.df)
  pc1.pc2 <- ggplot(pca.df,aes(x=PC1,y=PC2,color=group, fill=group)) +
    geom_point(size=4) +
    geom_polygon(data=hulls, aes(x=PC1, y=PC2), alpha=0.2)
  #pc1.pc2
  #ggsave(pc1.pc2,filename="pca.plot.PC1.PC2.png",width=6, height=4)


  tmp.df <- pca.df
  tmp.df$PC2 <- tmp.df$PC3
  hulls <- find_hulls(tmp.df)
  pc1.pc3 <- ggplot(pca.df,aes(x=PC1,y=PC3,color=group, fill=group)) +
    geom_point(size=4) +
    geom_polygon(data=hulls, aes(x=PC1, y=PC2), alpha=0.2)
  #pc1.pc3
  #ggsave(pc1.pc3,filename="pca.plot.PC1.PC3.png",width=6, height=4)

  tmp.df <- pca.df
  tmp.df$PC1 <- tmp.df$PC3
  hulls <- find_hulls(tmp.df)
  pc2.pc3 <- ggplot(pca.df,aes(x=PC2,y=PC3,color=group, fill=group)) +
    geom_point(size=4) +
    geom_polygon(data=hulls, aes(x=PC2, y=PC1), alpha=0.2)
  #pc2.pc3
  #ggsave(pc2.pc3,filename="pca.plot.PC2.PC3.png",width=6, height=4)

  pca.list <- list(p1 = pc1.pc2, p2 = pc1.pc3, pc3 = pc2.pc3)
  return(pca.list)

}

#' Principal Component Analysis Heatmap
#'
#' Generates heatmap using the PCA values as values
#' @param nano The nano object
#' @param groups Group information for the columns of the count matrix
#' @param countCutoff Mean cut value used to cut off unexpressed genes. Defaults to 5
#' @param cluster_row Logical, TRUE if rows should be clustered
#' @param cluster_col Logical, TRUE if columns should be clustered
#' @keywords principal component analysis pca heatmap
#' @export
#' @examples
#' plotPCAheatmap()
plotPCAheatmap <- function(nano, groups, countCutoff = 5, pcs = c(1:10), cluster_row = F, cluster_col = F){
  pcs <- pcs + 1
  # Extract expressed counts
  counts <- nano$counts
  no.cols <- ncol(counts)
  counts <- counts[,4:no.cols]
  mean.counts <- apply(X = counts, MARGIN = 1, FUN = mean)
  counts <- counts[mean.counts >= countCutoff,]

  # Calculate principal components
  pca.df <- createPCAdf(counts, groups)

  # Import pheatmap library
  library(pheatmap)

  pca.df.vals <- pca.df[,pcs]
  annot <- data.frame(Group = pca.df$group)
  pheatmap(t(pca.df.vals), annotation = annot, cluster_rows = cluster_row, cluster_cols = cluster_col)

}
