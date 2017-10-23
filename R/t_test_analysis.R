#' Generate numeric group vectors
#'
#' Generate numeric vectors characterizing the group columns
#' @param contrast The contrast for which the columns should be extracted
#' @param groups The groups the samples belong to
#' @keywords group column vector
#' @export
#' @examples
#' generateGroupVectors()
generateGroupVectors <- function(groups, contrast){
  g1 <- c()
  g2 <- c()
  group.labels <- strsplit(contrast, split="-")[[1]]
  for (i in 1:length(groups)){
    if (groups[i] == group.labels[1]){
      g1 <- c(g1, i)
    }
    else if (groups[i] == group.labels[2]){
      g2 <- c(g2, i)
    }
  }
  out.list <- list(g1 = g1, g2 = g2)
  return(out.list)
}

#' Fold change for a gene
#'
#' Calculate fold change of a gene between conditions
#' @param g1 Column vector of group 1
#' @param g2 Colum vector of group 2
#' @param gene Integer specifying the row of the gene
#' @param counts The count matrix
#' @keywords Fold change
#' @export
#' @examples
#' geneLogFC()
geneLogFC <- function(g1, g2, gene, counts){
  g1.mean <- mean(as.numeric(counts[gene, g1]))
  g2.mean <- mean(as.numeric(counts[gene, g2]))
  fc <- log2(g1.mean/g2.mean)
  return(fc)
}

#' Fold change for all genes
#'
#' Calculate fold change for all genes between a condition
#' @param counts The count matrix
#' @param contrast The contrast of interest
#' @param groups Character vectors specifying the groups of the columns
#' @keywords Fold change
#' @export
#' @examples
#' calcFCs()
calcFCs <- function(counts, contrast, groups){
  group.vecs <- generateGroupVectors(groups, contrast)
  g1 <- group.vecs$g1
  g2 <- group.vecs$g2
  
  fcs <- c()
  for(i in 1:nrow(counts)){
    fc <- geneLogFC(g1, g2, i, counts)
    fcs <- c(fcs, fc)
  }
  return(fcs)
}

#' Student t-Test
#'
#' Calculate P-value and t-Statistic using the t-Test
#' @param counts The count matrix
#' @param contrast The contrast of interest
#' @param groups Character vectors specifying the groups of the columns
#' @keywords Pvalue t-Test
#' @export
#' @examples
#' calcStudentPvals()
calcStudentPvals <- function(counts, contrast, groups){
  group.vecs <- generateGroupVectors(groups, contrast)
  g1 <- group.vecs$g1
  g2 <- group.vecs$g2
  
  pvals <- c()
  tstats <- c()
  for(i in 1:nrow(counts)){
    vec1 <- as.numeric(counts[i, g1])
    vec2 <- as.numeric(counts[i, g2])
    t <- t.test(vec1, vec2)
    pvals <- c(pvals, t$p.value)
    tstats <- c(tstats, t$statistic)
  }
  df <- data.frame(Pval = pvals, tStatistic = tstats)
  return(df)
}

#' Student t-Test Analysis
#'
#' Calculate fold changes, p-values and t-statistic for a contrast
#' @param counts The count matrix
#' @param contrast The contrast of interest
#' @param groups Character vectors specifying the groups of the columns
#' @keywords Pvalue t-Test t-Statistic fold change
#' @export
#' @examples
#' performTAnalysis()
performTAnalysis <- function(counts, contrast, groups, padj.method="BH"){
  logFCs <- calcFCs(counts, contrast, groups)
  stat.result <- calcStudentPvals(counts, contrast, groups)
  genes <- rownames(counts)
  # Adjust p values
  adj.p.vals <- p.adjust(as.numeric(stat.result$Pval), method=padj.method)
  
  df <- data.frame(Gene = genes, 
                   logFC = logFCs, 
                   P.Value = stat.result$Pval, 
                   adj.P.Val = adj.p.vals, 
                   tStatistic = stat.result$tStatistic)
  rownames(df) <- df$Gene
  return(df)
}
