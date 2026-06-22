# Internal: NA-fill km5 result structure used when clustering cannot proceed.
.makeNAKm5 <- function(n_rows, k) {
  bspart <- rep(NA_real_, n_rows)
  jac    <- rep(NA_real_, k)
  centr  <- rep(NA_real_, k)
  list(
    bspart    = bspart, jac  = jac, centr = centr,
    means     = bspart, bspart.or = bspart, bspart.inv = bspart,
    jac.or    = jac,    jac.inv   = jac,
    partition = bspart, jac.stab  = jac,
    csv = list(
      cluster_partition = NULL, cluster_mean         = NULL,
      cluster_centers   = NULL, cluster_size         = NULL,
      cluster_betweenss = NULL, cluster_totss        = NULL,
      cluster_tot.withinss = NULL, cluster_anova     = NULL
    )
  )
}

# Internal: process one metric across all k values for stability analysis.
# Called by both the serial (lapply) and parallel (parLapply) dispatchers.
# Returns list(stab_valores, csv_data).
.runStabilityForOneMetric <- function(i.metr, datos.bruto, names.metr, cbi, seed,
                                      k.range, bs, all_metrics, gold_standard,
                                      num.metrics, ...) {
  stab_valores <- rep(NA_real_, max(k.range))
  csv_data     <- vector("list", max(k.range))
  i            <- i.metr + 1L

  if (all_metrics) {
    message("Processing all metrics, 'merge', in dataframe (", length(names.metr), ")")
    data_to_cluster <- datos.bruto[, -1]
  } else {
    message("Processing metric: ", names.metr[i.metr], " (", i.metr, ")")
    data_to_cluster <- datos.bruto[, i]
  }

  for (j.k in k.range) {
    message("\tCalculation of k = ", j.k)
    v.size <- if (!all_metrics) length(levels(as.factor(data_to_cluster))) else j.k

    if (v.size >= j.k) {
      raw <- tryCatch(
        clusterbootWrapper(data = data_to_cluster, B = bs, bootmethod = "boot",
                           cbi = cbi, gold_standard = gold_standard,
                           krange = j.k, seed = seed, ...),
        error = function(e) e
      )
      if (inherits(raw, "error")) {
        message(paste0("\t", raw))
        message("\tWarning: Could not process data for k = ", j.k)
        km5 <- .makeNAKm5(length(datos.bruto[, i]), j.k)
      } else {
        km5         <- list()
        km5$cluster <- raw
        km5$jac     <- raw$bootmean
        km5$bspart  <- raw$partition
        km5$csv     <- list(
          cluster_partition    = raw$partition,
          cluster_mean         = raw$bootmean,
          cluster_centers      = raw$result$result$centers,
          cluster_size         = raw$result$result$size,
          cluster_betweenss    = raw$result$result$betweenss,
          cluster_totss        = raw$result$result$totss,
          cluster_tot.withinss = raw$result$result$tot.withinss,
          cluster_anova        = fAnova(raw$result$result, j.k, num.metrics)
        )
        km5$centr <- raw$result$result$centers
        for (km5.i in seq_along(km5$centr)) {
          km5$means[which(km5$bspart == km5.i)] <- km5$centr[km5.i]
        }
      }
    } else {
      message("\tWarning: Could not process data for k = ", j.k)
      km5 <- .makeNAKm5(length(datos.bruto[, i]), j.k)
    }

    stab_valores[j.k] <- mean(km5$jac)
    csv_data[[j.k]]   <- km5
  }

  list(stab_valores = stab_valores, csv_data = csv_data)
}

#' @title Stability index.
#' @name stability
#' @aliases stability
#' @description
#' This analysis permits to estimate whether the clustering is meaningfully
#' affected by small variations in the sample. First, a clustering using the
#' k-means algorithm is carried out. The value of \code{k} can be provided by the user.
#' Then, the stability index is the mean of the Jaccard coefficient
#' values of a number of \code{bs} bootstrap replicates. The values are in the range [0,1],
#' having the following meaning:
#' \itemize{
#' \item Unstable: [0, 0.60[.
#' \item Doubtful: [0.60, 0.75].
#' \item Stable: ]0.75, 0.85].
#' \item Highly Stable: ]0.85, 1].
#' }
#'
#' @param data A \code{\link{SummarizedExperiment}}.
#' The SummarizedExperiment must contain an assay with the following structure:
#' A valid header with names. The first  column of the header is the ID or name
#' of the instance of the dataset (e.g., ontology, pathway, etc.) on which the
#' metrics are measured.
#' The other columns of the header contains the names of the metrics.
#' The rows contains the measurements of the metrics for each instance in the dataset.
#' @param k Positive integer. Number of clusters between [2,15] range.
#' @param bs Positive integer. Bootstrap value to perform the resampling.
#' @param cbi Clusterboot interface name (default: "kmeans"):
#' "kmeans", "clara", "clara_pam", "hclust", "pamk", "pamk_pam", "pamk".
#' Any CBI appended with '_pam' makes use of \code{\link{pam}}.
#' The method used in 'hclust' CBI is "ward.D2".
#' @param getImages Boolean. If true, a plot is displayed.
#' @param all_metrics Boolean. If true, clustering is performed upon all the dataset.
#' @param seed Positive integer. A seed for internal bootstrap.
#' @param gold_standard Numeric vector. A vector of clusters from a gold standard classification, e.g. c(1,2,1,1,2).
#' Only applicable if parameter 'all_metrics' is set to TRUE.
#' @param numCores Number of cores to be used (>1 will use parallel processing)
#' @param ... Additional arguments passed to internal clustering functions.
#'
#' @return A \code{\link{ExperimentList}} containing the stability and cluster measurements
#'  for k clusters.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' result <- stability(ontMetrics, k=6, getImages=TRUE)
#'
#' @references
#' \insertRef{milligan1996measuring}{evaluomeR}
#'
#' \insertRef{jaccard1901distribution}{evaluomeR}
#'
#'
stability <- function(data, k=5, bs=100, cbi="kmeans",
                      getImages=FALSE, all_metrics=FALSE, seed=NULL,
                      gold_standard=NULL, numCores=1,...) {

  data <- assayAsDF(data)

  checkKValue(k)
  
  if(numCores>1){
    # Parallel version
    runStabilityIndex_parallel(data, k.min=k, k.max=k, bs=bs, cbi=cbi,
                      all_metrics=all_metrics, gold_standard=gold_standard, seed=seed, numCores=numCores, ...)
    stabilityDataFrame <- suppressWarnings(
      runStabilityIndexTableRange(data, k.min=k, k.max=k))
  }else{
    runStabilityIndex(data, k.min=k, k.max=k, bs=bs, cbi=cbi,
                      all_metrics=all_metrics, gold_standard=gold_standard, seed=seed, ...)
    stabilityDataFrame <- suppressWarnings(
      runStabilityIndexTableRange(data, k.min=k, k.max=k))
  }
  
  if (getImages == TRUE) {
    suppressWarnings(
      runStabilityIndexK_IMG(bs, k.min = k, k.max = k))
  }
  se <- createSEList(stabilityDataFrame)
  return(se)
}


#' @title Stability index for a range of k clusters.
#' @name stabilityRange
#' @aliases stabilityRange
#' @description
#' This analysis permits to estimate whether the clustering is meaningfully
#' affected by small variations in the sample. For a range of k values (\code{k.range}),
#' a clustering using the k-means algorithm is carried out.
#' Then, the stability index is the mean of the Jaccard coefficient
#' values of a number of \code{bs} bootstrap replicates. The values are in the range [0,1],
#' having the following meaning:
#' \itemize{
#' \item Unstable: [0, 0.60[.
#' \item Doubtful: [0.60, 0.75].
#' \item Stable: ]0.75, 0.85].
#' \item Highly Stable: ]0.85, 1].
#' }
#'
#' @inheritParams stability
#' @param k.range Concatenation of two positive integers.
#' The first value \code{k.range[1]} is considered as the lower bound of the range,
#' whilst the second one, \code{k.range[2]}, as the higher. Both values must be
#' contained in [2,15] range.
#'
#' @return A \code{\link{ExperimentList}} containing the stability and cluster measurements
#'  for 2 to \code{k} clusters.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' result <- stabilityRange(ontMetrics, k.range=c(2,3))
#'
#' @references
#' \insertRef{milligan1996measuring}{evaluomeR}
#'
#' \insertRef{jaccard1901distribution}{evaluomeR}
#'
#'
stabilityRange <- function(data, k.range=c(2,15), bs=100, cbi="kmeans",
                           getImages=FALSE, all_metrics=FALSE, seed=NULL,
                           gold_standard=NULL, numCores=1, ...) {
  k.range.length = length(k.range)
  if (k.range.length != 2) {
    stop("k.range length must be 2")
  }
  k.min = k.range[1]
  k.max = k.range[2]
  checkKValue(k.min)
  checkKValue(k.max)
  if (k.max < k.min) {
    stop("The first value of k.range cannot be greater than its second value")
  }
  data <- assayAsDF(data)
  
  if(numCores>1){
    # Parallel version
    runStabilityIndex_parallel(data, k.min=k.min, k.max=k.max, bs, cbi, all_metrics=all_metrics, gold_standard=gold_standard, seed=seed, numCores=numCores, ...)
    stabilityDataFrame <- suppressWarnings(
      runStabilityIndexTableRange(data, k.min=k.min, k.max=k.max))
  }else{
    runStabilityIndex(data, k.min=k.min, k.max=k.max, bs, cbi, all_metrics=all_metrics, gold_standard=gold_standard, seed=seed, ...)
    stabilityDataFrame <- suppressWarnings(
      runStabilityIndexTableRange(data, k.min=k.min, k.max=k.max))
  }

  if (getImages == TRUE) {
    suppressWarnings(
      runStabilityIndexK_IMG(bs, k.min=k.min, k.max=k.max))
    suppressWarnings(
      runStabilityIndexMetric_IMG(bs, k.min=k.min, k.max=k.max))
  }
  se <- createSEList(stabilityDataFrame)
  return(se)
}

#' @title Stability index for a set of k clusters.
#' @name stabilitySet
#' @aliases stabilitySet
#' @description
#' This analysis permits to estimate whether the clustering is meaningfully
#' affected by small variations in the sample. For a set of k values (\code{k.set}),
#' a clustering using the k-means algorithm is carried out.
#' Then, the stability index is the mean of the Jaccard coefficient
#' values of a number of \code{bs} bootstrap replicates. The values are in the range [0,1],
#' having the following meaning:
#' \itemize{
#' \item Unstable: [0, 0.60[.
#' \item Doubtful: [0.60, 0.75].
#' \item Stable: ]0.75, 0.85].
#' \item Highly Stable: ]0.85, 1].
#' }
#'
#' @inheritParams stability
#' @param k.set A list of integer values of \code{k}, as in c(2,4,8).
#' The values must be contained in [2,15] range.
#'
#' @return A \code{\link{ExperimentList}} containing the stability and cluster measurements
#'  of the list of \code{k} clusters.
#'
#' @examples
#' # Using example data from our package
#' data("rnaMetrics")
#' result <- stabilitySet(rnaMetrics, k.set=c(2,3))
#'
#' @references
#' \insertRef{milligan1996measuring}{evaluomeR}
#'
#' \insertRef{jaccard1901distribution}{evaluomeR}
#'
#'
stabilitySet <- function(data, k.set=c(2,3), bs=100, cbi="kmeans",
                         getImages=FALSE, all_metrics=FALSE, seed=NULL,
                         gold_standard=NULL, numCores=1, ...) {
  k.set.length = length(k.set)
  if (k.set.length == 0) {
    stop("k.set list is empty")
  } else if (k.set.length == 1) {
    stop("k.set list contains only one element. For one K analysis use 'stability' method")
  }
  k.set = sort(k.set)
  for (k in k.set) {
    checkKValue(k)
  }

  data <- assayAsDF(data)
  
  if(numCores>1){
    # Parallel version
    runStabilityIndex_parallel(data, k.set = k.set, bs=bs, cbi=cbi, all_metrics=all_metrics, gold_standard=gold_standard, seed=seed, numCores=numCores, ...)
    stabilityDataFrame <- suppressWarnings(
      runStabilityIndexTableRange(data, k.set = k.set))
  }else{
    runStabilityIndex(data, k.set = k.set, bs=bs, cbi=cbi, all_metrics=all_metrics, gold_standard=gold_standard, seed=seed, ...)
    stabilityDataFrame <- suppressWarnings(
      runStabilityIndexTableRange(data, k.set = k.set))
  }
  

  if (getImages == TRUE) {
    suppressWarnings(
      runStabilityIndexK_IMG(bs, k.set = k.set))
    suppressWarnings(
      runStabilityIndexMetric_IMG(bs, k.set = k.set))
  }
  se <- createSEList(stabilityDataFrame)
  return(se)
}

runStabilityIndex <- function(data, k.min=NULL, k.max=NULL, bs,
                              cbi, all_metrics, seed, k.set=NULL,
                              gold_standard=NULL,...) {
  if (is.null(seed)) {
    seed = pkg.env$seed
  }
  if (is.null(k.min) && is.null(k.max) && is.null(k.set)) {
    stop("runStabilityIndex: All k parameters are null!")
  }
  if (!is.null(gold_standard)) {
    if (all_metrics==FALSE) {
      stop("Gold standard parameter can be set only if the clustering of all the metrics is selected (all_metrics = TRUE)")
    }
    if (bs != 0) {
      message("Warning: 'gold_standard' parameter is set, argument 'bs' will be ignored.")
    }
  }

  data <- removeNAValues(data)
  dfStats(data)

  datos.bruto  <- data
  names.metr   <- names(datos.bruto)[-1]
  pkg.env$names.metr <- names.metr

  i.min <- k.min; i.max <- k.max
  if (!is.null(k.set)) {
    k.range <- k.set; k.range.length <- length(k.set)
  } else {
    k.range <- i.min:i.max; k.range.length <- length(i.min:i.max) + 1
  }

  num.metrics <- if (all_metrics) 1L else length(names.metr)
  if (all_metrics) pkg.env$names.metr <- "all_metrics"

  # --- serial dispatch via lapply, same body as parallel ---
  resultados <- lapply(seq_len(num.metrics), function(i.metr)
    .runStabilityForOneMetric(i.metr, datos.bruto, names.metr, cbi, seed,
                              k.range, bs, all_metrics, gold_standard, num.metrics, ...))

  m.stab.global     <- lapply(resultados, `[[`, "stab_valores")
  m.stab.global.csv <- lapply(resultados, `[[`, "csv_data")

  e.stab.global <- NULL
  for (j.k in k.range) {
    e.stab.global[[j.k]] <- matrix(data=NA, nrow=k.range.length, ncol=1)
    for (i.metr in seq_len(length(pkg.env$names.metr))) {
      e.stab.global[[j.k]][i.metr] <- m.stab.global[[i.metr]][j.k]
    }
  }

  pkg.env$m.stab.global.csv <- m.stab.global.csv
  pkg.env$m.stab.global     <- m.stab.global
  pkg.env$e.stab.global     <- e.stab.global
}

# Parallel version of runStabilityIndex
# Parallelized function that divides the metrics across workers (adds numCores)
runStabilityIndex_parallel <- function(data, k.min=NULL, k.max=NULL, bs,
                                       cbi, all_metrics, seed, k.set=NULL,
                                       gold_standard=NULL, numCores=NULL,...) {

  if (is.null(seed)) {
    seed = pkg.env$seed
  }
  if (is.null(k.min) && is.null(k.max) && is.null(k.set)) {
    stop("runStabilityIndex: All k parameters are null!")
  }
  if (!is.null(gold_standard)) {
    if (all_metrics==FALSE) {
      stop("Gold standard parameter can be set only if the clustering of all the metrics is selected (all_metrics = TRUE)")
    }
    if (bs != 0) {
      message("Warning: 'gold_standard' parameter is set, argument 'bs' will be ignored.")
    }
  }

  data <- removeNAValues(data)
  dfStats(data)

  datos.bruto  <- data
  names.metr   <- names(datos.bruto)[-1]
  pkg.env$names.metr <- names.metr

  i.min <- k.min; i.max <- k.max
  if (!is.null(k.set)) {
    k.range <- k.set; k.range.length <- length(k.set)
  } else {
    k.range <- i.min:i.max; k.range.length <- length(i.min:i.max) + 1
  }

  num.metrics <- if (all_metrics) 1L else length(names.metr)
  if (all_metrics) pkg.env$names.metr <- "all_metrics"

  message("Number of cores in parallelization: ", numCores)
  cl <- makeCluster(numCores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterEvalQ(cl, library(evaluomeR))

  # --- parallel dispatch: same helper as serial, different dispatcher ---
  resultados <- parLapply(cl, seq_len(num.metrics), function(i.metr)
    .runStabilityForOneMetric(i.metr, datos.bruto, names.metr, cbi, seed,
                              k.range, bs, all_metrics, gold_standard, num.metrics, ...))

  m.stab.global     <- lapply(resultados, `[[`, "stab_valores")
  m.stab.global.csv <- lapply(resultados, `[[`, "csv_data")

  e.stab.global <- NULL
  for (j.k in k.range) {
    e.stab.global[[j.k]] <- matrix(data=NA, nrow=k.range.length, ncol=1)
    for (i.metr in seq_len(num.metrics)) {
      e.stab.global[[j.k]][i.metr] <- m.stab.global[[i.metr]][j.k]
    }
  }

  pkg.env$m.stab.global.csv <- m.stab.global.csv
  pkg.env$m.stab.global     <- m.stab.global
  pkg.env$e.stab.global     <- e.stab.global
}

runStabilityIndexTableRange <- function(data, k.min=NULL, k.max=NULL, k.set=NULL) {
  if (is.null(k.min) && is.null(k.max) && is.null(k.set)) {
    stop("runStabilityIndexTableRange: All k parameters are null!")
  }
  stabilityDataFrame = NULL

  m.stab.global = pkg.env$m.stab.global
  m.stab.global.csv = pkg.env$m.stab.global.csv
  names.metr = pkg.env$names.metr

  measures = NULL
  # Key = Dataframe name - Value = Header name for each k
  measures["stability_mean"]= c("Mean_stability_k_")
  measures["cluster_partition"]= c("Cluster_partition_k_")
  measures["cluster_mean"]= c("Cluster_mean_k_")
  measures["cluster_centers"]= c("Cluster_centers_k_")
  measures["cluster_size"]= c("Cluster_size_k_")
  measures["cluster_betweenss"]= c("Cluster_betweenss_k_")
  measures["cluster_totss"]= c("Cluster_totss_k_")
  measures["cluster_tot.withinss"]= c("Cluster_tot.withinss_k_")
  measures["cluster_anova"]= c("Cluster_anova_k_")

  k.range = NULL
  k.range.length = NULL
  if (!is.null(k.set)) {
    k.range = k.set
    k.range.length = length(k.set)
  } else {
    k.range = k.min:k.max
    k.range.length = length(k.min:k.max)+1
  }

  for (measure in names(measures)) {
    stabilityDataList = list()
    # Build header
    header <- list("Metric")
    for (k in k.range) {
      nextHeader <- paste(measures[measure], k, sep="")
      header <- c(header, nextHeader)
    }
    header = unlist(header, use.names=FALSE)

    # Build rows
    for (i.metr in 1:length(names.metr)) {
      measure.data.list = list()
      wrapper=NULL # Wrapper object to extract the measure data
      wrapper$metric=names.metr[i.metr]
      if (measure == "stability_mean") {
        jacMean = m.stab.global[[i.metr]]
        cur.value = jacMean[c(k.range)]
        measure.data.list = c(measure.data.list, cur.value)
      } else {
        km5.list = m.stab.global.csv[[i.metr]]
        for (k in k.range) {
          km5.cur = km5.list[[k]]
          cur.value = getMeasureValue(km5.cur, measure)
          measure.data.list = c(measure.data.list, cur.value)
        }
      }
      wrapper[[measure]] = unlist(measure.data.list, use.names = FALSE)
      stabilityDataList[[i.metr]] = unlist(wrapper, use.names=FALSE)
    }  # end for i.metr

    # Transform into dataframe
    # Add df if has data
    if (ncol(t(data.frame(stabilityDataList))) != 1) {
      stabilityDataFrame[[measure]] = t(data.frame(stabilityDataList))
      colnames(stabilityDataFrame[[measure]]) = header
      rownames(stabilityDataFrame[[measure]]) <- NULL
    }
  }

  return(stabilityDataFrame)
}

# Stability index per K value (x values = matrics)
runStabilityIndexK_IMG <- function(bs, k.min = NULL, k.max = NULL, k.set = NULL) {
  if (is.null(k.min) && is.null(k.max) && is.null(k.set)) {
    stop("runStabilityIndexK_IMG: All k parameters are null!")
  }
  ancho=6
  alto=4
  escala=0.6
  escalax=escala
  ajuste=0.5
  escalat=0.5
  escalap=0.4
  contadorFiguras=1
  #Pattern: GlobalStability_K_2, ..., GlobalStability_K_N
  figurename="GlobalStability_K_"

  colores <- c("black", "blue", "red", "darkgreen", "orange", "magenta", "gray")
  ltype <- c(1, 2, 3, 4, 3, 4, 6)
  i.min=k.min
  i.max=k.max

  k.range = NULL
  k.range.length = NULL
  g.main = NULL

  if (!is.null(k.set)) {
    k.range = k.set
    k.range.length = length(k.set)
    setAsStrList = paste(as.character(k.set),collapse=", ",sep="")
    g.main=paste(" St. Indices of the metrics for k in {", setAsStrList, "}",sep="")
  } else {
    k.range = i.min:i.max
    k.range.length = length(k.min:k.max)+1
    g.main=paste(" St. Indices of the metrics for k in [", i.min, ",", i.max, "]",sep="")
  }

  stype <- c(1:k.range.length)

  e.stab.global = pkg.env$e.stab.global
  names.metr = pkg.env$names.metr
  margins <- par(mar=c(5,5,3,3))
  on.exit(par(margins))
  xnames=as.character(names.metr)
  ynames="Global Stability Indices"

  metrics_length = length(names.metr)

  num_metrics_plot = 19
  num_iterations = round(metrics_length/num_metrics_plot)
  if (num_iterations > 0) {
    num_iterations = num_iterations - 1
  }
  for (iteration in 0:num_iterations) {
    i = 1
    labels = list()
    rangeStart = (iteration*num_metrics_plot)+1
    rangeEnd = rangeStart+num_metrics_plot
    if (rangeEnd > metrics_length) {
      rangeEnd = metrics_length
    }
    new_xnames = xnames[rangeStart:rangeEnd]
    for (j.k in k.range) {
      cur.data = e.stab.global[[j.k]]
      cur.data = cur.data[rangeStart:rangeEnd]
      #cur.data = cur.data[!is.na(cur.data)]
      #cur.data = cur.data[c(i.min:i.max)
      plot(cur.data, main=g.main, axes=TRUE, col.axis="white",
           xlim=c(0.75,length(new_xnames)+0.25), xlab="", ylim=c(0,1),
           ylab=ynames, col="black", type="o", lwd=1, lty=stype[i], pch=stype[i])
      labels = c(labels,(paste0("k=", j.k)))
      i = i + 1
      par(new=TRUE)
    }
    axis(1,at=1:length(new_xnames),labels=new_xnames,las=2,cex.axis=0.75)
    axis(2,las=3,cex.axis=0.85)
    legend("bottomright", legend=labels, inset=.01, lwd=1, lty=stype,
           col="black", cex=0.7, pch=stype)
    mtext(side=1, text="Metrics",line=4)
    text(0.75, 0.9, "H.stab", cex=0.6, col = "black")
    abline(h = 0.85, col="black", lwd=1, lty=1) # Highly Stable: (0.85, 1]
    text(0.74, 0.8, "Stab", cex=0.6, col = "black")
    abline(h = 0.75, col="black", lwd=1, lty=1) # Stable: (0.75, 0.85]
    text(0.76, 0.65, "Doubt.", cex=0.6, col = "black")
    abline(h = 0.60, col="black", lwd=1, lty=1) # Doubtful: [0.60, 0.75]
    text(0.76, 0.05, "Unstab.", cex=0.6, col = "black")
    #abline(h = 0, col="black", lwd=1, lty=4) # Unstable: [0, 0.60)
    par(new=FALSE)
  }
}

# Stability index per metric (x values = k range)
runStabilityIndexMetric_IMG <- function(bs, k.min=NULL, k.max=NULL, k.set=NULL) {
  if (is.null(k.min) && is.null(k.max) && is.null(k.set)) {
    stop("runStabilityIndexMetric_IMG: All k parameters are null!")
  }
  ancho=6
  alto=4
  escala=0.9
  escalax=escala
  ajuste=0.5
  escalat=0.5
  escalap=0.4
  contadorFiguras=1
  #Pattern: GlobalStability_MetricX, ..., GlobalStability_MetricN
  figurename="GlobalStability_"

  colores <- c("black", "blue", "red", "darkgreen", "orange", "magenta", "gray")
  ltype <- c(1, 2, 3, 4, 3, 4, 6)

  stype <- c("o", "l")
  pchtype <- c(1, 2, 3, 4, 5, 5)


  k.range = NULL
  k.range.length = NULL
  if (!is.null(k.set)) {
    k.range = k.set
    k.range.length = length(k.set)
  } else {
    k.range = k.min:k.max
    k.range.length = length(k.min:k.max)
  }

  m.stab.global = pkg.env$m.stab.global
  names.metr = pkg.env$names.metr
  margins <- par(mar=c(5,5,3,3))
  on.exit(par(margins))
  for (i.metr in 1:length(names.metr)) {
    cur.data = m.stab.global[[i.metr]]
    cur.data = cur.data[!is.na(cur.data)]
    ymin = min(cur.data)
    if (is.na(ymin) || ymin == Inf) {
      ymin = 0
    }
    if (!is.null(k.set)) {
      setAsStrList = paste(as.character(k.set),collapse=", ",sep="")
      g.main=paste(" St. Indices of '", names.metr[i.metr], "' for k in {",
                   setAsStrList,"}",sep="")
    } else {
      g.main=paste(" St. Indices of '", names.metr[i.metr], "' for k in [",
                   k.min, ",", k.max,"]",sep="")
    }
    xnames=c(k.range)
    ynames="Global Stability Indices"

    plot(cur.data, main=g.main, axes=TRUE, col.axis="white",
         xlim=c(0.75,k.range.length+0.25), xlab="", ylim=c(ymin,1),
         ylab=ynames, col=colores[1],type="o", lwd=1, lty=ltype[1])
    axis(1,at=1:k.range.length,labels=xnames,las=1,cex.axis=escalax)
    axis(2,las=3,cex.axis=0.85)
    labels <- paste("b=", bs, sep = "")
    legend("topright", legend=labels, inset=.01, lwd=1, lty=ltype[1:4], col=colores[1:4], cex=0.7, pch=pchtype[1:4])
    mtext(side=1, text="K values",line=3)
  }

}