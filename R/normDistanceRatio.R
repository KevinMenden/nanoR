#' Distance Ratio Plot for all different normalizations
#'
#' Generates plot of distance ratio before and after normalization for each group
#' @param nano The non-normalized nano object
#' @param groups The groups the samples belong to
#' @keywords distance ratio normalization
#' @export
#' @examples
#' plotNormDistanceRatio()
plotNormDistanceRatio = function(nano, groups, countCutoff = 5, pcm = "mean", bm = "geo.mean"){
  # Background and positive control normalization
  nano <- nsBackgroundCorrect(nano, bm = bm)
  nano <- nsPositiveControlNormalization(nano, pcm = pcm)
  # Try out diferent normalizations
  nano.top100 <- nsNormalize(nano, "top100")
  nano.hk <- nsNormalize(nano, "housekeeping")
  nano.global <- nsNormalize(nano, "total")
  # Extract raw and normalized data
  rdata <- nano$bg.corr.counts
  top100.data <- nano.top100$counts
  hk.data <- nano.hk$counts
  global.data <- nano.global$counts
  
  # Remove annoation
  no.cols <- length(colnames(rdata))
  rdata <- removeAnnot(rdata)
  top100.data <- removeAnnot(top100.data)
  hk.data <- removeAnnot(hk.data)
  global.data <- removeAnnot(global.data)
  # Remove unexpressed
  rdata <- removeNonExpressed(rdata)
  top100.data <- removeNonExpressed(top100.data)
  hk.data <- removeNonExpressed(hk.data)
  global.data <- removeNonExpressed(global.data)
  
  
  # Attach group names to colnames
  colnames(rdata) <- paste(groups, colnames(rdata), sep="")
  colnames(top100.data) <- paste(groups, colnames(top100.data), sep="")
  colnames(hk.data) <- paste(groups, colnames(hk.data), sep="")
  colnames(global.data) <- paste(groups, colnames(global.data), sep="")

  
  # Calculate distance ratios
  raw.distance.ratio <- distanceRatio(rdata, groups)
  top100.distance.ratio <- distanceRatio(top100.data, groups)
  hk.distance.ratio <- distanceRatio(hk.data, groups)
  global.distance.ratio <- distanceRatio(global.data, groups)
  
  # Merge and melt dataframes
  levs <- levels(factor(groups))
  df.raw <- data.frame(groups = levs, values= raw.distance.ratio, variable=rep("Raw", length(levs)))
  df.top100 <- data.frame(groups = levs, values= top100.distance.ratio, variable=rep("Top100", length(levs)))
  df.hk <- data.frame(groups = levs, values= hk.distance.ratio, variable=rep("HK", length(levs)))
  df.global <- data.frame(groups = levs, values= global.distance.ratio, variable=rep("Global", length(levs)))
  mdata <- melt(rbind(df.raw, df.top100, df.hk, df.global))
  
  # Plot the Result
  p <- ggplot(mdata, aes(x=groups, y=value, fill=variable)) +
        geom_bar(stat="identity", position=position_dodge()) +
        labs( x = "Group", y = "Group to Global distance ratio", title = "Distance Ratio")
  
  return(p)
}
