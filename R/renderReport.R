rcc.dir <- "/sandata/11/NGSData/Transfer_OG/MasterKevin/Analysis_nanostring/"
design.file <-  "/sandata/11/NGSData/Transfer_OG/MasterKevin/Analysis_nanostring/design_ns_hdm.txt"

#' HTML Report Generation
#'
#' Function to generate a HTML report using Rmarkdown
#' @param rcc.dir The directory containing the RCC files
#' @param design.file The path to the design file
#' @param norm.method The method to use for RNA content normalization. Choose one of c("housekeeping", "top100", "total")
#' @param pcm Positive control method, the method to use for positive control normalization. Choose one of c("mean", "geo.mean")
#' @param bm Background method, the method to calculate the background values. Choose one of c("mean", "geo.mean")
#' @keywords HTML report
#' @export
#' @examples
#' renderReport()
renderReport <- function(rcc.dir, design.file, norm.method = "housekeeping", pcm = "geo.mean", bm = "mean"){
  library(rmarkdown)
  # Parse RCC files
  nano <- parseRCC(dir = rcc.dir)
  # Read design file
  design <- read.table(design.file, sep="\t", check.names=F, header=T)
  # Render the report
  # TODO: FIX PROBLEM WITH RMD FILE
  render('analysis_report.Rmd',
         params = list(nano = nano,
                       design = design,
                       norm.method = norm.method,
                       pcm = pcm,
                       bm = bm))

}


