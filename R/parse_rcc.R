#########################
## RCC file parser ######
#########################

# Note: The majority of the parsing code is a modified version of the RCC parsing function
# from the R Package NanoStringNorm


#' Parse single RCC file
#'
#' Function to parse the content of a single RCC file
#' @param path The file path
#' @param noGenes The number of genes in the panel
#' @keywords RCC parsing
#' @export
#' @examples
#' readSingleRCC()
readSingleRCC = function(path, noGenes){
  # Read the header, which is always 15 rows
  rcc.header  <- read.table(path, nrows=15, comment.char="<", sep=",", as.is=TRUE)
  # Read the actual count data
  rcc.data <- read.table(path, skip=25, header=TRUE, comment.char="<", as.is=TRUE,
                         sep=",", nrows=noGenes, check.names=FALSE)
  # Change "Count" column to $SAMPLE column
  sample <- basename(path)
  colnames(rcc.data)[4] <- sample
  colnames(rcc.header)[2] <- sample
  x <- list(rcc.data = rcc.data, rcc.header = rcc.header)
  return(x)
}

#' Parse RCC files
#'
#' Function to parse all RCC files from one directory
#' and merge them into a nano object
#' @param dir The path to the RCC directory. Defaults to the current directory.
#' @param rcc.pattern Pattern to recognize RCC files in case the don't end with .RCC, which is assumed by default.
#' @param walkSubdirs If TRUE, RCC files in subdirectories will be parsed as well. Default is FALSE.
#' @keywords RCC parsing
#' @export
#' @examples
#' parseRCC()
parseRCC = function(dir=".", rcc.pattern=".RCC$", walkSubdirs=FALSE){

  # List all RCC files in this directory and subdirectories
  rcc.files <- c()
  # Current directory
  new.files <- list.files(path=dir, pattern = rcc.pattern, full.names=TRUE)
  new.files <- sapply(new.files, normalizePath)
  rcc.files <- c(rcc.files, new.files)
  # Subdirectories
  if (walkSubdirs){
    for (subdir in list.dirs(recursive=FALSE)){
      new.files <- list.files(path=subdir, pattern = rcc.pattern)
      new.files <- sapply(new.files, normalizePath)
      rcc.files <- c(rcc.files, new.files)
    }
  }
  # Stop if no RCC files could be found
  if (length(rcc.files) == 0){
    stop("Could not find any RCC files!")
  }

  #testing
  lines <- readLines(rcc.files[1])
  lines <- lines[26:length(lines)]
  noGenes = 0
  tmp = 0
  for (l in lines){
    tmp = tmp + 1
    if (l == "</Code_Summary>"){
      noGenes = tmp
    }
  }
  noGenes <- noGenes - 3
  # Read and merge RCC files
  merged.headers <- NULL
  merged.data <- NULL
  count <- 1
  for (f in rcc.files){
    rcc <- readSingleRCC(f, noGenes)
    header <- rcc$rcc.header
    data <- rcc$rcc.data

    # First File
    if (count==1){
      merged.headers <- header
      merged.data <- data
    } else {
      # Rest of files
      merged.headers <- data.frame(merged.headers, subset(header, select = 2), check.names = F)
      merged.data <- data.frame(merged.data, subset(data, select = 4), check.names = F)
    }
    count = count + 1

  }
  # Some data adjustments
  merged.headers[3, 1] <- "sample.id"
  merged.headers[9, 1] <- "lane.id"
  colnames(merged.data[1]) <- "CodeCount"
  # Return data
  rownames(merged.data) <- merged.data$Name
  x <- list(counts = merged.data, header = merged.headers, raw.counts = merged.data)
  return(x)
}
