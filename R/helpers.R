#' @title Sample input data loader
#' @name loadSample
#' @aliases loadSample
#' @description
#' This method is a wrapper to load sample input data located inside evaluomeR package.
#'
#' @param descriptor Sample file to load: "ont-metrics", "rna-metrics" or "biopathways-metrics".
#'
#' @return The dataset specified via \code{descriptor} as a dataframe.
#'
#' @examples
#' # Using example data from our package
#' metrics = loadSample("rna-metrics")
#'
loadSample <- function(descriptor) {
  samples <- c('ont-metrics','rna-metrics','biopathways-metrics')
  if (is.element(descriptor, samples)) {
    dataFrame <- read.csv(file=system.file('extdata',descriptor, package="evaluomeR"), header=TRUE);
    return(dataFrame)
  } else {
    stop("Invalid descriptor")
  }
}

#' @title Dataframe getter for \code{qualityRange} function.
#' @name getDataQualityRange
#' @aliases getDataQualityRange
#' @description
#' This method is a wrapper to retrieve a specific dataframe given a \code{k} value from
#' the object returned by \code{\link{qualityRange}} function.
#'
#' @param data The object returned by \code{\link{qualityRange}} function.
#' @param k The desired \code{k} cluster.
#'
#' @return The dataframe that contains information about the selected \code{k} cluster.
#'
#' @examples
#' # Using example data from our package
#' metrics = loadSample("ont-metrics")
#' qualityRangeData <- qualityRange(data=metrics, k.range=c(3,5), getImages = FALSE)
#' # Getting dataframe that contains information about k=5
#' k5DataFrame = getDataQualityRange(qualityRangeData, 5)
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

#' @title SummarizedExperiment to Dataframe
#' @name seToDataFrame
#' @aliases seToDataFrame
#' @description
#' This method is a wrapper to transform a SummarizedExperiment object to a
#' Dataframe processable in our methods.
#'
#' @param SummarizedExperiment A \code{SummarizedExperiment} object
#' (see \code{\link{SummarizedExperiment}}).
#'
#' @return The dataframe that contains information of the first
#' assay in \code{SummarizedExperiment}.
#'
#' @examples
#' # Using example data from airway package
#' library(airway)
#' data(airway)
#' airwayData = seToDataFrame(airway)
#' airwayData = airwayData[1:10000,1:4]
#' stability(airwayData, bs = 20, getImages=FALSE)
#' correlations(airwayData, getImages=FALSE)
#'
seToDataFrame <- function(SummarizedExperiment) {
  se=SummarizedExperiment
  if (length(assays(se)) == 0) {
    stop("SummarizedExperiment has no assays, length is 0")
  }
  test = assay(se,1)
  Datasets <- NULL
  if (is.null(rownames(test))) {
    Datasets <- paste("Dataset_", c(1:length(test[,1])), sep="")
  } else {
    Datasets <- rownames(test)
  }

  test <- data.frame(Datasets,test)
  return(test)
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
