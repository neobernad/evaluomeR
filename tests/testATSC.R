library(evaluomeR)
library(RSKC)

data("golub")
dataset = golub

# First, data cleaning

r_cleanDataset = evaluomeR::cleanDataset(dataset, correlation_threshold = 1)
dataset =  r_cleanDataset$dataset
R =  r_cleanDataset$R


pca_suitability = evaluomeR::PCASuitability(R, sig_level = 0.05)
print(pca_suitability$pca_suitable)

if (pca_suitability$pca_suitable) {
  message("PCA is suitable")
  r_pca = evaluomeR::performPCA(dataset)
  dataset = r_pca$dataset_ncp
}

head(dataset)

# Second clustering and optimal k

ATSC <- function(data, k.range=c(2,15), bs=100, cbi="kmeans",
                 a_max = 0.1, all_metrics=TRUE, seed=NULL) {
  k.range.length = length(k.range)
  if (k.range.length != 2) {
    stop("k.range length must be 2")
  }
  k.min = k.range[1]
  k.max = k.range[2]
  #checkKValue(k.min)
  #checkKValue(k.max)
  if (k.max < k.min) {
    stop("The first value of k.range cannot be greater than its second value")
  }
  data = as.data.frame(SummarizedExperiment::assay(data))
  stabRange = evaluomeR::stabilityRange(data=data, cbi=cbi, k=k.range, bs=bs,
                             all_metrics = all_metrics, seed=seed)
  stab = evaluomeR::standardizeStabilityData(stabRange, k.range = k.range)
  qualRange = evaluomeR::qualityRange(data=data, cbi=cbi, k=k.range,
                             all_metrics = all_metrics, seed = seed)
  qual = evaluomeR::standardizeQualityData(qualRange, k.range = k.range)
  rOptimalK = evaluomeR::getOptimalKValue(stabRange, qualRange)
  optimalK = rOptimalK$Global_optimal_k
  message(paste0("Optimal k: ", optimalK))

  # Automated Trimmed & Sparse Clustering
  L1 = evaluomeR::getRSKCL1Boundry(data, k=optimalK, seed=seed)
  alpha = evaluomeR::getRSKCAlpha(data, k=optimalK, L1=L1, seed=seed)

  print(message(paste0("L1: ", L1)))
  print(message(paste0("alpha: ", alpha)))
  return (list(stab=stab, qual=qual))

}



test = ATSC(data=dataset, k.range=c(2,6), cbi="clara")

#stab <- stabilityRange(data=test, cbi="rskc", alpha=0, L1=2, k=k.range, bs=100, all_metrics = TRUE)

ncol(dataset[-1])

rskc_out = RSKC(dataset[-1], L1=8, alpha=0, ncl=3)
#union_vector = c(rskc_out$oE,rskc_out$oW)
#union_vector_unique = unique(union_vector)
#union_vector_unique = sort(union_vector_unique)

names_w_zero = names(rskc_out$weights)[rskc_out$weights == 0]
dataset_cleaned <- dataset[, !(names(dataset) %in% names_w_zero)]
#union_vector = c(test$oE,test$oW)
#print(union_vector)
#length(unique(union_vector))

#dataset_cleaned <- dataset[-union_vector, ]

