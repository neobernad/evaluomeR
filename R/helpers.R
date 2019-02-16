#' @title Sample input data loader
#' @name loadSample
#' @aliases loadSample
#' @description
#' This method is a wrapper to load sample input data located inside evaluomeR package.
#'
#' @param descriptor Sample file to load: "ont-metrics", "rna-metrics" or "biopathways-metrics".
#'
#' @return The \code{\link{SummarizedExperiment}} specified via \code{descriptor}.
#'
#' @examples
#' # Using example data from our package
#' metrics = loadSample("rna-metrics")
#'
loadSample <- function(descriptor) {
  samples <- c('ont-metrics','rna-metrics','biopathways-metrics')
  if (is.element(descriptor, samples)) {
    dataFrame <- read.csv(file=system.file('extdata',descriptor, package="evaluomeR"), header=TRUE);
    se <- createSE(dataFrame)
    return(se)
  } else {
    stop("Invalid descriptor")
  }
}

#' @title Dataframe getter for \code{qualityRange} function.
#' @name getDataQualityRange
#' @aliases getDataQualityRange
#' @description
#' This method is a wrapper to retrieve a specific \code{\link{SummarizedExperiment}} given a \code{k} value from
#' the object returned by \code{\link{qualityRange}} function.
#'
#' @param data The object returned by \code{\link{qualityRange}} function.
#' @param k The desired \code{k} cluster.
#'
#' @return The \code{\link{SummarizedExperiment}} that contains information about the selected \code{k} cluster.
#'
#' @examples
#' # Using example data from our package
#' metrics = loadSample("ont-metrics")
#' qualityRangeData <- qualityRange(data=metrics, k.range=c(3,5), getImages = FALSE)
#' # Getting dataframe that contains information about k=5
#' k5Data = getDataQualityRange(qualityRangeData, 5)
#'
getDataQualityRange <- function(data, k) {
  dataNames = names(data)
  # From, e.g, k_4 to 4
  kValues = substr(dataNames, nchar(dataNames)-0, nchar(dataNames))
  kValues = as.integer(kValues)
  kValues.length = length(kValues);

  if (k >= kValues[1] && k <= kValues[kValues.length]) {
    column = paste("k_", k, sep="")

    return(data[[column]])
  } else {
    error=paste("Selected k (",k,") is not in the range of k.range ["
                , kValues[1], ",", kValues[kValues.length], "]", sep="")
    stop(error)
  }
}

#####################
## Private methods ##
#####################

# Source: http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
# Supresses the output messages of a function 'x'
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

isString <- function(str) {
  if (is.character(str) & length(str) > 0) {
    return(TRUE)
  }
  error=paste("Input parameter '", str, "' is not a string", sep="")
  stop(error)
}

## Inserts a newrow to existingDF at position r
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

## Checks if k is in range [2,15]
checkKValue <- function(k) {
  if (k < 2 || k > 15) {
    stop("k value is not in range [2,15]")
  }
}

checkDirectory <- function(path) {
  if (is.null(path) || !dir.exists(path)) {
    error <- paste("Invalid path '", path, "'. Directory does not exists or it is not a directory", sep="")
    stop(error)
  }
  if (!endsWith(path, "/")) {
    path <- paste(path, "/", sep="")
  }
  return(path)
}

getAssay <- function(SummarizedExperiment, position) {
  se=SummarizedExperiment
  se.length <- length(assays(se))
  if (se.length == 0) {
    stop("SummarizedExperiment has no assays, length is 0")
  }
  if (position > se.length) {
    error <- paste("SummarizedExperiment has no assay in position ",
                   position, sep="")
    stop(error)
  }
  test = assay(se, position)
  # Datasets <- test[,1]
  # if (is.null(rownames(test))) {
  #   Datasets <- paste("Dataset_", c(1:length(test[,1])), sep="")
  # } else {
  #   Datasets <- rownames(test)
  # }

  # test <- data.frame(Datasets,test)
  test <- data.frame(test)
  names(test) <- colnames(test)
  return(test)
}

# data: One dataframe, thus one assay
createSE <- function(data) {
  nrows <- nrow(data); ncols <- ncol(data)
  counts <- data.matrix(data)
  colnames(counts) <- NULL
  colData <- DataFrame(metrics=colnames(data),
                       row.names=colnames(data))
  se <- SummarizedExperiment(assays=SimpleList(counts),
                              colData=colData)
  return(se)
}

# data: A list of dataframes
createSEList <- function(data) {
  if (!is.list(data)) {
    stop("Input variable is not a list")
  }
  if (length(data) == 0) {
    stop("Input variable is an empty list")
  }
  length = length(names(data))
  seList <- list()
  for (i in 1:length) {
    cur.data <- data[[i]]
    dataMatrix <- suppressWarnings(data.matrix(cur.data))
    dataMatrix[,1] <- cur.data$Metric
    se <- createSE(dataMatrix)
    seList <- c(seList, se)
  }
  names(seList) <- names(data)
  return(seList)
}
