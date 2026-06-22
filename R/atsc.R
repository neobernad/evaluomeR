ATSC <- function(data, k.range=c(2,15), bs=100, cbi="clara",
                 max_alpha = 0.1,
                 L1=NULL, alpha=NULL, gold_standard=NULL,
                 seed=NULL, numCores=1, clusteringSparsity="kmeans") {
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
  if (is.null(seed)) {
    seed = pkg.env$seed
  }
  all_metrics=TRUE
  message(paste0("Computing optimal k value with '", cbi, "'"))

  data = assayAsDF(data)
  # Stability indexes
  stabRange = evaluomeR::stabilityRange(data=data, cbi=cbi, k=k.range, bs=bs,
                                        all_metrics = all_metrics,
                                        gold_standard=gold_standard, seed=seed, numCores=1)
  stab = evaluomeR::standardizeStabilityData(stabRange, k.range = k.range)
  # Quality indexes
  qualRange = evaluomeR::qualityRange(data=data, cbi=cbi, k=k.range,
                                      all_metrics = all_metrics, seed = seed, numCores=1)
  # Optimal k analysis
  qual = evaluomeR::standardizeQualityData(qualRange, k.range = k.range)
  rOptimalK = evaluomeR::getOptimalKValue(stabRange, qualRange)
  optimalK = as.numeric(rOptimalK$Global_optimal_k)
  message(paste0("Optimal k: ", optimalK))

  # Automated Trimmed & Sparse Clustering
  message("Determining best L1 and alpha parameter automatically, it might take a while...")
  if (is.null(L1)) {
    L1 = evaluomeR::getRSKCL1Boundry(data, k=optimalK, clustering=clusteringSparsity, seed=seed)
  }

  if (is.null(alpha)) {
      best_alphas = evaluomeR::getRSKCAlpha(data, k=optimalK, L1=L1,
                                    max_alpha = max_alpha, seed=seed, numCores=numCores)
      alpha = best_alphas$best_alpha_qual
  }
  #invisible(suppressMessages({}))


  message(paste0("\tUsing L1 '", L1, "' and alpha '", alpha, "'"))

  message("Running Trimmed & Sparse Clustering algorithm")
  set.seed(seed)
  rskcOut = RSKC(data[, !colnames(data) %in% "Description"], L1=L1, alpha=alpha, ncl=optimalK)
  data_trimmed = data
  trimmedRows = c()
  gold_standard_trimmed = c()
  if (alpha > 0 && length(rskcOut$oW) > 1) { # Rows were trimmed, detect which
    trimmedRows  = c(rskcOut$oE,rskcOut$oW)
    trimmedRows = unique(trimmedRows)
    trimmedRows = sort(trimmedRows)
    message(paste0("\tNumber of trimmed cases: ", length(trimmedRows)))
    data_trimmed = data_trimmed[-trimmedRows, ]
    if (!is.null(gold_standard)) { # Remove trimmed cases from gold standard vector
      message("\tTrimming gold standard as well")
      gold_standard = gold_standard[-trimmedRows]
      gold_standard_trimmed = gold_standard
    }

  }
  # Check if L1 removed columns (feature weight == 0)
  trimmedColumns = names(rskcOut$weights)[rskcOut$weights == 0]
  if (length(trimmedColumns) > 0) {
    message(paste0("\tNumber of affected columns: ", length(trimmedColumns)))
    data_trimmed = data_trimmed[, !(names(data_trimmed) %in% trimmedColumns)]
    #message(paste(trimmedColumns, sep = ","))
  }

  message("Computing optimal k value on the dataset processed by a trimmed sparse clustering method.")
  # Repeat optimal k analysis for 'data_trimmed'
  stabRange_ATSC = evaluomeR::stabilityRange(data=data_trimmed, cbi=cbi, k=k.range, bs=bs,
                                             all_metrics = all_metrics,
                                             gold_standard=gold_standard, seed=seed, numCores=1)
  stab_ATSC = evaluomeR::standardizeStabilityData(stabRange_ATSC, k.range = k.range)

  qualRange_ATSC = evaluomeR::qualityRange(data=data_trimmed, cbi=cbi, k=k.range,
                                           all_metrics = all_metrics, seed = seed, numCores=1)
  qual_ATSC = evaluomeR::standardizeQualityData(qualRange_ATSC, k.range = k.range)

  rOptimalK_ATSC = evaluomeR::getOptimalKValue(stabRange_ATSC, qualRange_ATSC)
  optimalK_ATSC = as.numeric(rOptimalK_ATSC$Global_optimal_k)
  message(paste0("New optimal k with ATSC: ", optimalK_ATSC))



  return (list(
    # Before ATSC
    stab=stab, qual=qual, optimalK=optimalK,
    # After ATSC
    stab_ATSC=stab_ATSC, qual_ATSC=qual_ATSC, optimalK_ATSC=optimalK_ATSC,
    # Additional parameters of interes
    rskcOut=rskcOut, trimmedRows=trimmedRows, trimmedColumns=trimmedColumns,
    trimmmedDataset=data_trimmed, L1=L1, alpha=alpha,
    gold_standard_trimmed=gold_standard_trimmed
  ))

}
