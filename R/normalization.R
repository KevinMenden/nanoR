#############################
### Normalization Section ###
#############################
#' Normalization functions
#'
#' Content normalization of nanostring data
#' @param nsdata A nano object
#' @param method The normalization method. One of c("top100", "housekeeping", "total")
#' @keywords normalization
#' @export
#' @examples
#' nsNormalize()
nsNormalize = function(nsdata, method){
  # Get counts and number of samples
  counts <- nsdata$counts
  no.cols <- length(colnames(counts))
  lane.norm.factors <- c()
  endo.counts <- counts[grepl("Endogenous", counts$CodeClass), 4:no.cols]
  names <- counts$Name[grepl("Endogenous", counts$CodeClass)]

  # Top 100 normalization
  if (method == "top100"){
    average.expression <- apply(
      X = endo.counts,
      MARGIN = 1,
      FUN = mean
    )
    ae.genes <- data.frame(gene = names, ae = average.expression)
    ae.genes <- ae.genes[with(ae.genes, order(ae, decreasing = TRUE)),]
    top100 <- as.character(ae.genes$gene[1:100])
    top100.counts <- counts[counts$Name %in% top100,4:no.cols]
    mean.top100 <- apply(
      X = top100.counts,
      MARGIN = 2,
      FUN = geo.mean
    )
    global.mean.top100 <- mean(mean.top100)
    lane.norm.factors <- global.mean.top100 / mean.top100
  }

  # Housekeeping normalization
  if (method == "housekeeping"){
    hk.counts <- counts[counts$CodeClass == "Housekeeping",4:no.cols]
    mean.hk <- apply(
      X = hk.counts,
      MARGIN = 2,
      FUN = geo.mean
    )
    global.mean.hk <- mean(mean.hk)
    lane.norm.factors <- global.mean.hk / mean.hk
  }

  # Total sum
  if (method == "total"){
    total.counts <- endo.counts
    sum.total <- apply(
      X = total.counts,
      MARGIN = 2,
      FUN = sum
    )
    total.mean <- mean(sum.total)
    lane.norm.factors <- total.mean/sum.total
  }

  # Low CV methods (experimental)
  if (method == "lowcv"){
    mean.genes <- apply(
      X = endo.counts,
      MARGIN = 1,
      FUN = mean
    )
    endo.counts.at <- endo.counts[mean.genes >= 5,]
    cv.list <- calculateMatrixCV(endo.counts.at)
    norm.counts <- endo.counts.at[order(cv.list),]
    norm.counts <- norm.counts[1:50,]

    mean.lowcv <- apply(
      X = norm.counts,
      MARGIN = 2,
      FUN = geo.mean
    )
    global.mean.lowcv <- mean(mean.lowcv)
    lane.norm.factors <- global.mean.lowcv / mean.lowcv

  }

  # multiply with factors
  x <- t(apply(
    X = counts[,4:no.cols],
    MARGIN = 1,
    FUN = '*',
    lane.norm.factors
  ))
  counts <- cbind(counts[,1:3],x)
  nsdata$counts <- counts
  nsdata$norm.factors <- lane.norm.factors
  return(nsdata)
}



#' Positive Control Normalization
#'
#' Normalize using the positive controls (spike-ins)
#' @param nsdata A nano object
#' @param pcm Positive control method, should be one of c("mean", "geo.mean")
#' @keywords positive control normalization
#' @export
#' @examples
#' nsPositiveControlNormalization()
nsPositiveControlNormalization <- function(nsdata, pcm = "geo.mean"){
  counts <- nsdata$counts
  no.cols <- length(colnames(counts))
  pos.counts <- counts[counts$CodeClass == "Positive",4:no.cols]
  pos.means.samples <- apply(
    X = pos.counts,
    MARGIN = 2,
    FUN = pcm
  )
  pos.mean <- mean(pos.means.samples)
  pos.factors <- pos.mean / pos.means.samples
  # multiply with factors
  x <- t(apply(
    X = counts[,4:no.cols],
    MARGIN = 1,
    FUN = '*',
    pos.factors
  ))
  counts <- cbind(counts[,1:3],x)
  nsdata$counts <- counts
  nsdata$pos.factors <- pos.factors
  return(nsdata)
}

#' Background Correction
#'
#' Calculates background noise and subtracts it from counts
#' @param nsdata A nano object
#' @param sd.factor Number of standard deviations to be added to background value. Defaults to 2
#' @param bm Background method, currently only "mean" possible
#' @keywords background correction negative controls
#' @export
#' @examples
#' nsBackgroundCorrect()
# Currently: use the mean + 2sd as background and subtract that from every count value
nsBackgroundCorrect = function(nsdata, bm = "mean", sd.factor = 2){
  counts <- nsdata$counts
  negative.counts <- counts[counts$CodeClass == "Negative", ]
  no.cols <- length(colnames(counts))
  bg_values <- c()
  # Calculate background values for each lane
  for (i in 4:no.cols){
    lane <- negative.counts[,i]
    neg.sd <- sd(lane)
    neg.mean <- mean(lane)
    neg.value <- neg.mean + (sd.factor*neg.sd)
    bg_values <- c(bg_values, neg.value)
  }

  # Subtract background value from each count value
  x <- t(apply(
    X = counts[,4:no.cols],
    FUN = '-',
    MARGIN = 1,
    bg_values
  ))
  x[x < 0] <- 1
  nsdata$raw.counts <- nsdata$counts
  nsdata$counts <- cbind(counts[,1:3],x)
  nsdata$bg.corr.counts <- nsdata$counts
  nsdata$background <- bg_values
  return(nsdata)
}


