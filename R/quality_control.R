#################################
### Quality Control Functions ###
#################################

#' Binding Density Plot
#'
#' Plot the binding densities for every sample for QC
#' @param nsdata A nano object
#' @keywords binding density qc
#' @export
#' @examples
#' plotBindingDensities()
plotBindingDensities = function(nsdata){
  header <- nsdata$header
  sample.names <- colnames(header)
  sample.names <- sample.names[2:length(sample.names)]
  binding.densities <- c()
  bd_cutoff_low <- 0.05
  bd_cutoff_high <- 2.25
  for (i in 2:length(colnames(header))){
    bd <- header[14,i]
    binding.densities <- c(binding.densities, bd)
  }
  binding.densities <- as.numeric(binding.densities)
  
  f = function(x){x >= bd_cutoff_low && x <= bd_cutoff_high}
  pass <- sapply(binding.densities, f)
  qc.pass <- rep("True", length(binding.densities))
  qc.pass[pass] <- "True"
  qc.pass[!pass] <- "False"
  
  df <- data.frame(samples = sample.names, binding.densities, Pass = qc.pass)
  
  p <- ggplot(data=df, aes(x=samples, y=binding.densities, fill=Pass)) +
    geom_bar(stat="identity", color="blue") +
    theme(axis.text.x=element_text(size=5, angle = 90, vjust = 0.5)) +
    xlab("Sample") +
    ylab("Binding Density") +
    geom_hline(yintercept = bd_cutoff_low, colour="purple") +
    geom_hline(yintercept = bd_cutoff_high, colour="purple") +
    scale_fill_manual(values = c(True = "lightblue", False = "orange"))
  p
}


#' Field of View count plot
#'
#' Plot the field of view count for every sample for QC
#' @param nsdata A nano object
#' @keywords fov qc
#' @export
#' @examples
#' plotFOV()
plotFOV = function(nsdata){
  header <- nsdata$header
  sample.names <- colnames(header)
  sample.names <- sample.names[2:length(sample.names)]
  fov_cutoff = 0.75
  fovs <- c()
  for (i in 2:length(colnames(header))){
    fov.all <- as.numeric(header[10,i])
    fov.actual <- as.numeric(header[11,i])
    fov.percent <- fov.actual / fov.all
    fovs <- c(fovs, fov.percent)
  }

  qc.pass <- rep("True", length(fovs))
  qc.pass[fovs >= fov_cutoff] <- "True"
  qc.pass[!fovs >= fov_cutoff] <-"False"
  
  df <- data.frame(sample.names, fovs, Pass = qc.pass)

  p <- ggplot(data=df, aes(x=sample.names, y=fovs, fill = Pass)) +
    geom_bar(stat="identity", color="blue") +
    theme(axis.text.x=element_text(size=5, angle = 90, vjust = 0.5)) +
    xlab("Sample") +
    ylab("Percentage FOV counted") +
    geom_hline(yintercept = fov_cutoff, colour="purple") +
    scale_fill_manual(values = c(True = "lightblue", False = "orange"))
  p
}

#' Positive scaling factor plot
#'
#' Plot the scaling factors calculated from the positive controls
#' @param nsdata A nano object
#' @keywords positive control factor
#' @export
#' @examples
#' plotPositiveScalingFactors()
plotPositiveScalingFactors = function(nsdata){
  counts <- nsdata$counts
  no.cols <- length(colnames(counts))
  counts <- counts[counts$CodeClass == "Positive",]
  ps_cutoff_low <- 0.3
  ps_cutoff_high <- 3.0
  spikes <- c()
  for (i in 4:no.cols){
    m <- mean(counts[,i])
    spikes <- c(spikes,m)
  }
  mean.spikes <- mean(spikes)
  scaling_factors <- c()
  for (i in 1:length(spikes)){
    sf <- mean.spikes/spikes[i]
    scaling_factors <- c(scaling_factors, sf)
  }
  samples <- colnames(counts)[4:no.cols]
  
  f = function(x){x >= ps_cutoff_low && x <= ps_cutoff_high}
  pass <- sapply(scaling_factors, f)
  qc.pass <- rep("True", length(scaling_factors))
  qc.pass[pass] <- "True"
  qc.pass[!pass] <- "False"
  
  df <- data.frame(samples, scaling_factors, Pass = qc.pass)
  
  p <- ggplot(data=df, aes(x=samples, y=scaling_factors, fill=Pass)) +
    geom_bar(stat="identity", color="blue") +
    theme(axis.text.x=element_text(size=5, angle = 90, vjust = 0.5)) +
    xlab("Sample") +
    ylab("Positive Scaling Factor") +
    geom_hline(yintercept = ps_cutoff_low, colour = "purple") +
    geom_hline(yintercept = ps_cutoff_high, colour = "purple") +
    scale_fill_manual(values = c(True = "lightblue", False = "orange"))
  p
}

#' Normalization factor plot
#'
#' Plot the calculated normalization factors
#' @param nsdata A nano object
#' @keywords normalization factor
#' @export
#' @examples
#' plotNormFactors()
plotNormFactors = function(nsdata){
  counts <- nsdata$counts
  scaling_factors <- nsdata$norm.factors
  samples <- colnames(counts)[4:ncol(counts)]
  nf_cutoff_low <- 0.1
  nf_cutoff_high <- 10.0
  
  f = function(x){x >= nf_cutoff_low && x <= nf_cutoff_high}
  pass <- sapply(scaling_factors, f)
  qc.pass <- rep("True", length(scaling_factors))
  qc.pass[pass] <- "True"
  qc.pass[!pass] <- "False"
  
  df <- data.frame(samples, scaling_factors, Pass = qc.pass)

  p <- ggplot(data=df, aes(x=samples, y=scaling_factors, fill=Pass)) +
    geom_bar(stat="identity", color="blue") +
    theme(axis.text.x=element_text(size=5, angle = 90, vjust = 0.5)) +
    xlab("Sample") +
    ylab("Content Normalization Factor") +
    geom_hline(yintercept = nf_cutoff_low, colour = "purple") +
    geom_hline(yintercept = nf_cutoff_high, colour = "purple") +
    scale_fill_manual(values = c(True = "lightblue", False = "orange"))
  p
}

#' Background value plot
#'
#' Plot the calculated background values
#' @param nsdata A nano object
#' @keywords background threshold
#' @export
#' @examples
#' plotBackground()
plotBackgrund = function(nsdata){
  
  bg.vals <- nsdata$background
  if (is.null(bg.vals)){
    nano2 <- nsdata
    nano2 <- nsBackgroundCorrect(nano2)
    bg.vals <- nano2$background
  }
  
  samples <- colnames(nsdata$counts)[4:ncol(nsdata$counts)]
  
  df <- data.frame(samples, bg.vals)
  
  p <- ggplot(data=df, aes(x=samples, y=bg.vals)) +
    geom_bar(stat="identity", color="blue", fill = "lightblue") +
    theme(axis.text.x=element_text(size=5, angle = 90, vjust = 0.5)) +
    xlab("Sample") +
    ylab("Background Threshold Value") +
    geom_hline(yintercept = nf_cutoff_low, colour = "purple") +
    geom_hline(yintercept = nf_cutoff_high, colour = "purple") 
  p
}




#' Housekeeping genes plot
#'
#' Plot the variation of housekeeping genes
#' @param nano The nano object
#' @keywords housekeeping linechart
#' @export
#' @examples
#' plotHousekeepingGenes()
plotHousekeepingGenes = function(nano){
  counts <- nano$counts
  hk.counts <- counts[counts$CodeClass == "Housekeeping",]
  mdata <- melt(hk.counts)

  ggplot(data = mdata, aes(x=variable, y=value, group=Name, shape=Name, colour=Name)) +
    geom_line(aes(linetype=Name), size=1) + geom_point(size=2, fill="white") +
    ggtitle("Housekeeping Genes Expression") +
    ylab("Count") +
    xlab("Sample") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

}
#' Correlation between housekeeping genes
#'
#' Plots a correlation matrix for the housekeeping genes.
#' @param nano The nano object
#' @keywords housekeeping correlation
#' @export
#' @examples
#' housekeepingCorrelation()
housekeepingCorrelation = function(nano){
  counts <- nano$counts
  hk.counts <- counts[counts$CodeClass == "Housekeeping",]
  hk.counts <- removeAnnot(hk.counts)
  hk.cor <- cor(t(hk.counts), method="pearson")
  mcormat <- melt(hk.cor)
  ggplot(data = mcormat, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradient(limit = c(0,1), space = "Lab",name="Pearson\nCorrelation")
}


#' Correlation Matrix Top 100 Genes
#'
#' Plot correlation matrix for the top 100 genes
#' @param nano The nano object
#' @keywords top100 correlation
#' @export
#' @examples
#' topGenesCorrelation()
topGenesCorrelation = function(nano){
  endo.counts <- nano$counts[grepl("Endogenous", nano$counts$CodeClass),]
  endo.counts <- removeAnnot(endo.counts)
  average.expression <- apply(
    X = endo.counts,
    MARGIN = 1,
    FUN = mean
  )
  endo.counts <- endo.counts[order(average.expression, decreasing = TRUE),]
  top100.counts <- endo.counts[1:100,]
  top.cor <- cor(t(top100.counts), method="pearson")
  mcormat <- melt(top.cor)
  ggplot(data = mcormat, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(limit = c(-1,1), space = "Lab",name="Pearson\nCorrelation") +
    theme(axis.text.x = element_blank())
}


#' Plot R2 of positive control counts
#'
#' Plot correlation matrix for the top 100 genes
#' @param nano The nano object
#' @keywords R2 positive control linear
#' @export
#' @examples
#' positiveControlR2()
positiveControlR2 = function(nano){
  r2.cutoff = 0.95
  counts <- nano$counts
  pos <- counts[counts$CodeClass == "Positive",]
  pos <- pos[,4:ncol(pos)]
  concs <- c(0.125, 0.5, 2, 8, 32, 128)
  # Put in increasing order
  order.vals <- apply(pos, 1, mean)
  pos <- pos[order(order.vals),]
  r2.f = function(x){summary(lm(x~concs))$r.squared}
  r2.vals <- apply(pos, 2, r2.f)
  
  
  qc.pass <- rep("True", length(r2.vals))
  qc.pass[r2.vals >= r2.cutoff] <- "True"
  qc.pass[!r2.vals >= r2.cutoff] <-"False"
  
  df <- data.frame(Sample = colnames(pos), R2 = r2.vals, Pass = qc.pass)
  
  p <- ggplot(data=df, aes(x=Sample, y=R2, fill=Pass)) +
    geom_bar(stat="identity", color="blue") +
    theme(axis.text.x=element_text(size=5, angle = 90, vjust = 0.5)) +
    xlab("Sample") +
    ylab("R2") +
    geom_hline(yintercept = r2.cutoff, colour = "purple") +
    scale_fill_manual(values = c(True = "lightblue", False = "orange"))
  p
  
}




