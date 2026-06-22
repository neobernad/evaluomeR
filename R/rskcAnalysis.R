annotateClustersByMetric <- function(df, k.range, bs, seed){
  if (is.null(seed)) {
    seed = pkg.env$seed
  }
  df <- assayAsDF(df)
  # Create a dataframe by removing NAs from the original data.
  df_clean = na.omit(df)

  # Compute stability, quality and the optimal k value from evaluome
  stabilityData <- stabilityRange(data=df_clean, k.range=k.range,
                                  bs=bs, getImages = FALSE, seed=seed)

  qualityData <- qualityRange(data=df_clean, k.range=k.range,
                              getImages = FALSE, seed=seed)

  kOptTable <- getOptimalKValue(stabilityData, qualityData)

  # Get the clusters obtained by evaluome for each k, together with
  # the optimal k
  clusters = as.data.frame(assay(stabilityData$cluster_partition))
  clusters$optimal_k = kOptTable$Global_optimal_k

  # Compose the results
  # For each metric, get the optimal k, get the clusters formed by using
  # that optimal k, include this information in a dataframe
  result_list = list()

  for (i in 1:nrow(clusters)){
    metric = as.character(clusters$Metric[i])
    # Get the optimal k
    optimal_k = clusters$optimal_k[i]

    # Get the clusters formed by using the optimal k as a vector of integers
    optimal_cluster = dplyr::select(clusters, dplyr::contains(as.character(optimal_k)))
    optimal_cluster = optimal_cluster[i,1]
    optimal_cluster = as.character(optimal_cluster)
    optimal_cluster = as.numeric(strsplit(optimal_cluster, ", ")[[1]])

    # Create a dataframe including the individual id, the concerning metric
    # and the cluster id in which the individual is classfied.
    annotated_df_clean = dplyr::select(df_clean, 1)
    annotated_df_clean$cluster = optimal_cluster

    # Merge this dataframe with the original one, so that original individuals
    # removed due to NAs will be present an NA as cluster.
    # Include this dataframe in the named list, using the name of the metric as
    # a key.
    result_list[[metric]] = merge(dplyr::select(df, 1, dplyr::contains(metric)), annotated_df_clean, all.x = TRUE)
  }
  return(result_list)
}


#' @title Get the range of each metric per cluster from the optimal cluster.
#' getMetricRangeByCluster
#' @aliases getMetricRangeByCluster
#' @description
#' Obtains the ranges of the metrics obtained by each optimal cluster.
#'
#' @param df Input data frame. The first column denotes the identifier of the
#' evaluated individuals. The remaining columns contain the metrics used to
#' evaluate the individuals. Rows with NA values will be ignored.
#' @param k.range Range of k values in which the optimal k will be searched
#' @param bs Bootstrap re-sample param.
#' @param seed Random seed to be used.
#'
#' @return A dataframe including the min and the max value for each
#' pair (metric, cluster).
#' @export
#'


getMetricRangeByCluster <- function(df, k.range, bs, seed) {
  if (is.null(seed)) {
    seed = pkg.env$seed
  }
  df <- assayAsDF(df)
  annotated_clusters_by_metric = annotateClustersByMetric(df, k.range, bs, seed)

  metrics = c()
  cluster_ids = c()
  min_values = c()
  max_values = c()
  # For each metric
  for (metric in names(annotated_clusters_by_metric)){
    # Get the optimal clusters
    annotated_clusters = annotated_clusters_by_metric[[metric]]

    # For each cluster, get the minimal and the maximal value
    for (cluster_id in 1:max(annotated_clusters$cluster, na.rm=T)) {
      concrete_cluster_values = dplyr::filter(annotated_clusters, cluster==cluster_id) %>% dplyr::pull(metric)
      metrics = c(metrics, metric)
      cluster_ids = c(cluster_ids, cluster_id)
      min_values = c(min_values, min(concrete_cluster_values))
      max_values = c(max_values, max(concrete_cluster_values))
    }
  }
  return(data.frame(metric=metrics, cluster=cluster_ids, min_value=min_values, max_value=max_values))
}


#' @title Get the relevancy of each metric.
#' @name getMetricsRelevancy
#' @aliases getMetricsRelevancy
#' @description
#' Obtains the relevancy of the metrics using RSKC.
#'
#' @param df Input data frame. The first column denotes the identifier of the
#' evaluated individuals. The remaining columns contain the metrics used to
#' evaluate the individuals. Rows with NA values will be ignored.
#' @param k K value (number of clusters)
#' @param alpha 0 <= alpha <= 1, the proportion of the cases to be trimmed in robust sparse K-means, see \code{\link{RSKC}}.
#' @param L1 A single L1 bound on weights (the feature weights), see \code{\link{RSKC}}.
#' @param seed Random seed to be used.
#'
#' @return A dataframe including the min and the max value for each
#' pair (metric, cluster).
#' @export
#'
#' @examples
#' data("ontMetrics")
#' metricsRelevancy = getMetricsRelevancy(ontMetrics, k=3, alpha=0.1, seed=100)
#' metricsRelevancy$rskc # RSKC output object
#' metricsRelevancy$trimmed_cases # Trimmed cases from input (row indexes)
#' metricsRelevancy$relevancy # Metrics relevancy table
#'
getMetricsRelevancy <- function(df, k, alpha=0, L1=NULL, seed=NULL) {
    if (is.null(seed)) {
    seed = pkg.env$seed
  }
#  if (is.null(alpha)) {
#    alpha = 0.1
#  }

  df <- assayAsDF(df)

  if (is.null(L1)) {
    print(paste0("No L1 provided. Computing best L1 boundry with 'sparcl::KMeansSparseCluster.permute'"))
    dataMatrix = as.matrix(df)
    wbounds = seq(2,sqrt(ncol(dataMatrix)), len=30)
    km.perm <- sparcl::KMeansSparseCluster.permute(dataMatrix,K=k,wbounds=wbounds,nperms=5,silent=TRUE)
    L1 = km.perm$bestw

  }


  # Compute RSKC
  print(paste0("Alpha set as: ", alpha))
  print(paste0("L1 set as: ", L1))
  rskc_out = RSKC(df, k, alpha, L1 = L1, nstart = 200,
                  silent=TRUE, scaling = FALSE, correlation = FALSE)
  # Get trimmed cases from input
  union_vector = c(rskc_out$oE,rskc_out$oW)
  union_vector_unique = unique(union_vector)
  union_vector_unique = sort(union_vector_unique)

  # Metrics relevancy
  columns = c('metric', 'weight')
  rskc_df = data.frame(matrix(ncol = length(columns), nrow = length(rskc_out$weights)))
  colnames(rskc_df) = columns
  rskc_df['metric'] = names(rskc_out$weights)
  rskc_df['weight'] = rskc_out$weights
  rskc_df_sorted = rskc_df[order(rskc_df$weight, decreasing = TRUE), ] # Sorting from greater values to lower



  output = NULL
  output$rskc = rskc_out
  output$trimmed_cases = union_vector_unique
  output$relevancy = rskc_df_sorted
  return (output)
}


#' @title Computes L1 boundry
#' getRSKCL1Boundry
#' @aliases getRSKCL1Boundry
#' @description
#' Computes the L1 boundry for an RSKC cbi execution.
#'
#' @param df Input data frame. The first column denotes the identifier of the
#' evaluated individuals. The remaining columns contain the metrics used to
#' evaluate the individuals. Rows with NA values will be ignored.
#' @param k K value (number of clusters), only used when clustering = "kmeans" 
#' @param clustering Clustering method to be used
#' @param seed Random seed to be used.
#'
#' @return A single L1 bound on weights (the feature weights), see \code{\link{RSKC}}.
#' @export
#'
#' @examples
#' data("ontMetrics")
#' K-means (default)
#' l1_boundry = getRSKCL1Boundry(ontMetrics, k=3, clustering="kmeans", seed=100, silent=TRUE)
#' Hierarchical
#' l1_boundry = getRSKCL1Boundry(ontMetrics, clustering="hierarchical", seed=100)
#' 
getRSKCL1Boundry <- function(df, k, clustering="kmeans", seed=NULL) {
  if (is.null(seed)) {
    seed = pkg.env$seed
  }

  df <- assayAsDF(df)
  df_data = df[-1] # Removing 'Description' column as it is not numeric
  dataMatrix = as.matrix(df_data)
  wbounds = seq(2,sqrt(ncol(dataMatrix)), len=30)
  
  if(clustering == "kmeans"){
    message(paste0("Computing best L1 boundry with 'sparcl::KMeansSparseCluster.permute'"))
    km.perm <- sparcl::KMeansSparseCluster.permute(dataMatrix,K=k,wbounds=wbounds,nperms=5,silent=TRUE)
    #L1 = floor(km.perm$bestw)
    
  } else if (clustering == "hierarchical"){
    message(paste0("Computing best L1 boundry with 'sparcl::HierarchicalSparseCluster.permute'"))
    km.perm <- sparcl::HierarchicalSparseCluster.permute(dataMatrix,wbounds=wbounds,nperms=5)
    #L1 = floor(km.perm$bestw)
  }
  
  L1 = km.perm$bestw
  message(paste0("Best L1 found is: ", L1))
  return (as.numeric(L1))
}

#' @title Computes best alpha using Robust and Sparse Clustering
#' getRSKCAlpha
#' @aliases getRSKCAlpha
#' @description
#' Computes the proportion of the cases to be trimmed in robust sparse K-means, 0 <= alpha <= 1, see \code{\link{RSKC}}.
#'
#' @param df Input data frame. The first column denotes the identifier of the
#' evaluated individuals. The remaining columns contain the metrics used to
#' evaluate the individuals. Rows with NA values will be ignored.
#' @param k K value (number of clusters)
#' @param L1 A single L1 bound on weights (the feature weights), see \code{\link{RSKC}}.
#' @param max_alpha Maximum value of alpha,  iterating over seq(0, max_alpha, 0.05). Default is 0.1.
#' @param seed Random seed to be used.
#' @param numCores Number of cores to be used (>1 will use parallel processing)
#'
#' @return Best suitable alpha.
#' @export
#'
#' @examples
#' data("ontMetrics")
#' alpha = getRSKCAlpha(ontMetrics, k=3, L1=2, seed=100)
#'
getRSKCAlpha <- function(df, k, L1, max_alpha = 0.1, seed=NULL, numCores=1) {
  if (is.null(seed)) {
    seed = pkg.env$seed
  }

  df <- assayAsDF(df)

  # Helper structures
  ## Structure to keep track of stability and qualit of each RSKC run.
  rskc_run <- function(stab, qual, alpha) {
    structure(list(stab = stab, mean_stab = mean(as.double(stab)),
                   qual = qual, mean_qual = mean(as.double(qual)),
                   alpha = alpha), class = "rskc_run")
  }
  run_list = list()
  ## Getting best RSKC run according to stabilities and qualities
  best_stab = -Inf
  best_qual = -Inf
  best_alpha_stab = 0
  best_alpha_qual = 0
  getBestRun <- function(run_list) {
    for (run in run_list) {
      if (run$mean_stab > best_stab) {
        best_stab = run$mean_stab
        best_alpha_stab = run$alpha
      }
      if (run$mean_qual > best_qual) {
        best_qual = run$mean_qual
        best_alpha_qual = run$alpha
      }
    }

    return(
      list("best_alpha_stab" = best_alpha_stab,
           "best_alpha_qual"=best_alpha_qual,
           "best_stab"=best_stab,
           "best_qual"=best_qual)
    )
  }
  
  alpha_values = seq(0, max_alpha, 0.01)
  num_alphas = length(alpha_values)
  
  if(numCores>1){
    # The maximum parallel processes equals the number of alpha values to evaluate
    if(numCores > num_alphas){
      numCores <- num_alphas;
    }
    
    message("Number of cores in parallelization: ", numCores)
    cl <- makeCluster(numCores)
    on.exit(stopCluster(cl), add = TRUE)
    # Export variables from the current environment to the cluster
    clusterExport(cl, c("df", "k", "L1", "seed"), envir=environment())
    
    # Cargamos el paquete de evaluomeR
    clusterEvalQ(cl, {
      library(evaluomeR)
      library(RSKC)
    })
    run_list <- parLapply(cl, alpha_values, function(alpha){
      message(paste0("Running stability and quality indexes with alpha=",
                     alpha," k=",k, " L1=", L1))
      invisible(suppressMessages({
        stab = stability(data=df, k=k,
                         bs=100, seed=seed,
                         all_metrics=TRUE,
                         cbi="rskc", L1=L1, alpha=alpha)
        stab_table = standardizeStabilityData(stab)
        
        qual = quality(data=df, k=k,
                       seed=seed,
                       all_metrics=TRUE,
                       cbi="rskc", L1=L1, alpha=alpha)
        qual_table = standardizeQualityData(qual)
      }))
      
      # Evitamos exportar funciones a los procesos paralelos para reducir overhead
      return(structure(list(stab = stab_table, mean_stab = mean(as.double(stab_table)),
                     qual = qual_table, mean_qual = mean(as.double(qual_table)),
                     alpha = alpha), class = "rskc_run"))
    })
  }else{
    
    index = 1
    for (alpha in alpha_values) {
      message(paste0("Running stability and quality indexes with alpha=",
                     alpha," k=",k, " L1=", L1))
      invisible(suppressMessages({
        stab = stability(data=df, k=k,
                         bs=100, seed=seed,
                         all_metrics=TRUE,
                         cbi="rskc", L1=L1, alpha=alpha)
        stab_table = standardizeStabilityData(stab)
  
        qual = quality(data=df, k=k,
                       seed=seed,
                       all_metrics=TRUE,
                       cbi="rskc", L1=L1, alpha=alpha)
        qual_table = standardizeQualityData(qual)
      }))
  
      run_list[[index]] <- rskc_run(stab = stab_table, qual = qual_table, alpha=alpha)
      index = index + 1
    }
  }

  best_run = getBestRun(run_list)
  message(paste0("Highest stability found when alpha=", best_run$best_alpha_stab," (", best_run$best_stab,")"))
  message(paste0("Highest quality found when alpha=", best_run$best_alpha_qual," (", best_run$best_qual,")"))

  best_alpha = best_run$best_alpha_qual
  if (best_run$best_alpha_stab <= best_run$best_alpha_qual) {
    best_alpha = best_run$best_alpha_stab
  }
  #message(paste0("Using alpha=", best_alpha, " as it trims less data."))

  return (list(best_alpha_stab=as.numeric(best_run$best_alpha_stab),
               best_alpha_qual=as.numeric(best_run$best_alpha_qual)))
}

#' @title Automated Trimmed & Sparse Clustering
#' @name ATSC
#' @aliases ATSC
#' @description Automated Trimmed & Sparse Clustering.
#' This methods performs an optimal k value analysis with \code{\link{stabilityRange}}, \code{\link{qualityRange}}
#' and \code{\link{getOptimalKValue}} evaluomeR methods.
#' The optimal \code{k} value is used to compute estimate a \code{L1} bound and an \code{alpha} trimming portion automatically
#' in order to perform an automatic trimmed and sparse clustering.
#' This posibily results in the input dataset being trimmed (either by columns, determined by \code{L1} or
#' by rows, determined by \code{alpha}).
#' Another optimal k value analysis is then executed over the trimmed dataset, to conclude with the an optimal partition.
#'
#'
#' @param data A \code{\link{SummarizedExperiment}}.
#' The SummarizedExperiment must contain an assay with the following structure:
#' A valid header with names. The first column of the header is the ID or name
#' of the instance of the dataset (e.g., ontology, pathway, etc.) on which the
#' metrics are measured.
#' The other columns of the header contains the names of the metrics.
#' The rows contains the measurements of the metrics for each instance in the dataset.
#' @param k.range Concatenation of two positive integers.
#' The first value \code{k.range[1]} is considered as the lower bound of the range,
#' whilst the second one, \code{k.range[2]}, as the higher. Both values must be
#' contained in [2,15] range.
#' @param bs Positive integer. Bootstrap value to perform the resampling.
#' @param cbi Clusterboot interface name (default: "clara"):
#' "kmeans", "clara", "clara_pam", "hclust", "pamk", "pamk_pam", "rskc".
#' @param max_alpha Maximum value of alpha, iterating over seq(0, max_alpha, 0.01).
#' @param L1 A single L1 bound on weights (the feature weights), see \code{\link{RSKC}}.
#' If NULL, it is computed automatically.
#' @param alpha Trimming portion for RSKC, see \code{\link{RSKC}}.
#' If NULL, it is computed automatically via \code{max_alpha}.
#' @param gold_standard Numeric vector. A vector of clusters from a gold standard
#' classification. Only used when \code{all_metrics = TRUE}.
#' @param seed Positive integer. A seed for internal bootstrap.
#' @param numCores Number of cores to be used for the alpha search (>1 will use parallel processing).
#' @param clusteringSparsity Clustering method used to determine the L1 bound (default: "kmeans").
#'
#' @return A list containing:
#' \item{stab}{A data frame containing standardized stability.}
#' \item{qual}{A data frame containing standardized quality.}
#' \item{optimalK}{The optimal k value representing the optimal number of clusters determined from the initial analysis.}
#' \item{stab_ATSC}{A data frame containing standardized stability after applying ATSC.}
#' \item{qual_ATSC}{A data frame containing standardized quality  applying ATSC.}
#' \item{optimalK_ATSC}{The optimal k value representing the optimal number of clusters determined after applying ATSC.}
#' \item{rskcOut}{An object returned by the RSKC function containing clustering results, including weights and trimmed observations.}
#' \item{trimmedRows}{A vector of indices representing the rows that were trimmed from the dataset during the clustering process.}
#' \item{trimmedColumns}{A vector of names representing the columns that were trimmed (i.e., removed) from the dataset due to zero weights.}
#' \item{trimmedDataset}{A data frame containing the final processed dataset after trimming rows and columns.}
#' @export
#'
