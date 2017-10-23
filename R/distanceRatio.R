#' Distance Ratio Calculation
#'
#' Calculates list of distance ratios for every group in a count matrix
#' @param mat The count matrix
#' @param groups The groups the samples belong to
#' @keywords distance ratio
#' @export
#' @examples
#' distanceRatio()
distanceRatio = function(mat, groups){
  
  groups <- levels(factor(groups))
  ratio_list = c()
  
  for(g in groups){
    # Mean group distance
    sub_mat <- mat[,grepl(g, colnames(mat))]
    d <- dist(t(sub_mat))
    mean_group_distance <- mean(d)
    # Mean global distance
    d.glob <- dist(t(mat))
    d.glob <- as.matrix(d.glob)
    d.glob <- d.glob[grepl(g, rownames(d.glob)),]
    d.glob <- d.glob[, !grepl(g, colnames(d.glob))]
    mean_global_distance <- mean(d.glob)
    
    # Calculate distance ratio
    distance_ratio <- mean_group_distance/mean_global_distance
    
    ratio_list = c(ratio_list, distance_ratio)
  }
  ratio_list
}


#' Distance Ratio Plot
#'
#' Generates plot of distance ratio before and after normalization for each group
#' @param mat The nano object
#' @param groups The groups the samples belong to
#' @keywords distance ratio
#' @export
#' @examples
#' plotDistanceRatio()
plotDistanceRatio = function(mat, groups){
  # Extract raw and normalized data
  rdata <- mat$bg.corr.counts
  ndata <- mat$counts
  no.cols <- length(colnames(rdata))
  rdata <- rdata[,4:no.cols]
  ndata <- ndata[,4:no.cols]
  mean.counts <- apply(X = ndata, MARGIN = 1, FUN = mean)
  rdata <- rdata[mean.counts >= 5,]
  ndata <- ndata[mean.counts >= 5,]
  
  # Attach group names to colnames
  cols <- colnames(rdata)
  new_cols <- paste(groups, cols, sep="")
  colnames(rdata) <- new_cols
  colnames(ndata) <- new_cols
  
  # Calculate distance ratios
  raw.distance.ratio <- distanceRatio(rdata, groups)
  norm.distance.ratio <- distanceRatio(ndata, groups)
  
  # Merge and melt dataframes
  levs <- levels(factor(groups))
  df.raw <- data.frame(groups = levs, values= raw.distance.ratio, variable=rep("Raw", length(levs)))
  df.norm <- data.frame(groups = levs, values= norm.distance.ratio, variable=rep("Norm", length(levs)))
  mdata <- melt(rbind(df.raw, df.norm))
  
  # Plot the Result
  ggplot(mdata, aes(x=groups, y=value, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge()) +
    xlab("Group") +
    ylab("Group to Global distance ratio")
  
}
