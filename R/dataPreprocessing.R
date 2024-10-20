#' @title Computes L1 boundry
#' getRSKCL1Boundry
#' @aliases getRSKCL1Boundry
#' @description
#' Computes the L1 boundry for an RSKC cbi execution.
#'
#' @param df Input data frame. The first column denotes the identifier of the
#' evaluated individuals. The remaining columns contain the metrics used to
#' evaluate the individuals. Rows with NA values will be ignored.
#' @param k K value (number of clusters)
#' @param seed Random seed to be used.
#'
#' @return A single L1 bound on weights (the feature weights), see \code{\link{RSKC}}.
#' @export
#'
#' @examples
#' data("ontMetrics")
#' l1_boundry = getRSKCL1Boundry(ontMetrics, k=3, seed=100)
#'
getRSKCL1Boundry <- function(df, k, seed=NULL) {
  if (is.null(seed)) {
    seed = pkg.env$seed
  }

  df <- as.data.frame(assay(df))
  df_data = df[-1] # Removing 'Description' column as it is not numeric

  message(paste0("Computing best L1 boundry with 'sparcl::KMeansSparseCluster.permute'"))
  dataMatrix = as.matrix(df_data)
  wbounds = seq(2,sqrt(ncol(dataMatrix)), len=30)
  km.perm <- sparcl::KMeansSparseCluster.permute(dataMatrix,K=k,wbounds=wbounds,nperms=5,silent=TRUE)
  L1 = floor(km.perm$bestw)
  message(paste0("Best L1 found is: ", km.perm$bestw, ", using floor: ", L1))

  return (L1)
}
