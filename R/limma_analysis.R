###########################
### Limma Analysis ########
###########################
library(limma)

#' limma differential expression analysis
#'
#' Performs differential expression analysis on the data using the limma package
#' @param nano The nano object
#' @param designfile A design file containin the columns file, group, include
#' @param generatePlots Logical whether plots should be generated as well
#' @param meanCountCutoff The mean count value necessary for a gene to be considered for analysis. Defaults to 5
#' @keywords limma dge differential expression
#' @export
#' @examples
#' calcDELimma()
calcDELimma = function(nano, designfile, generatePlots = TRUE, meanCountCutoff = 5){
  design <- read.table(designfile, sep="\t", comment.char="", header=T)


  # Order design file and counts matrix, so that they have the same ordering
  # Delete samples with include == 0 in design file
  design <- design[design$include == 1,]
  counts <- nano$counts[grepl("Endogenous", nano$counts$CodeClass),]
  counts <- removeAnnot(counts)
  counts <- counts[,colnames(counts) %in% design$file]
  design <- design[with(design, order(file)),]
  counts <- counts[, with(counts, order(colnames(counts)))]
  #nano$counts <- cbind(nano$counts[,1:3], counts)
  no.cols <- ncol(counts)
  counts.log2 <- log2(counts+1)

  #Define the groups and contrasts
  groups <- design$group
  conts <- design$comparison[design$comparison != ""]


  # Consider only genes with a decent average expression
  mean.counts <- apply(X = counts, MARGIN = 1, FUN = mean)
  expressed <- mean.counts >= meanCountCutoff
  counts.expressed <- counts[expressed,]
  logcounts.expressed <- counts.log2[expressed,]

  logcounts.matrix <- as.matrix(logcounts.expressed)
  counts.matrix <- as.matrix(counts.expressed)


  # Only generate plots if specified
  if (generatePlots){
    # Create distance ration plot
    dist.plot <- plotDistanceRatio(nano, groups)
    ggsave(dist.plot, filename="Distance_ratio_plot.pdf")

    # Create PCAs
    plotPCA(nano, groups)

    # Log Count Distribution Plot
    png("Log2_count_distribution.png")
    hist(logcounts.matrix, breaks=20, main="Log2 count distribution", xlab="Log2 Counts", col="lightgreen", freq=FALSE)
    curve(dnorm(x, mean=mean(logcounts.matrix), sd=sd(logcounts.matrix)), add=TRUE, col="darkblue", lwd=2)
    dev.off()
    # Count Distribution Plot
    png("Count_distribution.png")
    hist(counts.matrix, breaks=20, main="Count distribution", xlab="Counts", col="lightgreen", freq=FALSE)
    curve(dnorm(x, mean=mean(counts.matrix), sd=sd(counts.matrix)), add=TRUE, col="darkblue", lwd=2)
    dev.off()
  }


  # Perform test for normality
  if (length(as.numeric(unlist(logcounts.matrix))) <= 5000){
    print("SHAPIRO WILK TEST FOR NORMALITY")
    shapiro.test(as.numeric(unlist(logcounts.matrix)))
  } else {
    print("Kolmogorov-Smirnov test")
    ks.test(x=as.numeric(unlist(logcounts.matrix)), y='pnorm', alternative='two.sided')
  }

  # Linear model fitting and DE analysis
  G <- factor(groups)
  dm <- model.matrix(~ -1 + G)
  colnames(dm) <- levels(G)

  contrasts <- makeContrasts(contrasts = conts, levels = dm)

  fit <- lmFit(logcounts.expressed, dm)
  fit2 <- contrasts.fit(fit = fit, contrasts = contrasts)
  fit2 <- eBayes(fit2)
  res <- decideTests(fit2)
  summary(res)

  # Make volcano plots and
  # Write separate table for every comparison
  print("Generating volcano plots ...")
  # Transfer conts to character vector, since it can cause trouble as factor
  conts <- as.character(conts)
  for (cont in conts){
    res <- topTable(fit2, coef=cont, number=Inf, sort.by="P")
    name <- paste("DGE_", cont, ".txt", sep="")
    write.table(res, name, sep="\t", quote=F, col.names = NA)

    res$threshold <- as.factor(abs(res$logFC) > 1 & res$adj.P.Val < 0.05)
    res$comparison <- cont
    res$geneID <- rownames(res)

    volcanoPlot(dge.result = res, contrast = cont)
  }


}
