library(evaluomeR)
library(RSKC)

#data("ontMetricsOBO")
#dataset = ontMetricsOBO

data("golub")
dataset = golub

# Gold standard vector

gold_standard_class = as.factor(as.vector(unlist(golub["Class"])))
gold_standard_class

level_mapping <- c("B" = 1, "T" = 2, "M" = 3)
map_strings_to_numbers <- function(strings) {
  factorized <- factor(strings, levels = names(level_mapping))
  return (as.numeric(factorized))
}

gold_standard_vector = as.numeric(lapply(gold_standard_class, map_strings_to_numbers))


# First, data cleaning

r_cleanDataset = evaluomeR::cleanDataset(dataset, correlation_threshold = 1)
dataset =  r_cleanDataset$dataset

pca_suitability = evaluomeR::PCASuitability(r_cleanDataset$R, sig_level = 0.05)

if (pca_suitability$pca_suitable) {
  message("PCA is suitable")
  r_pca = evaluomeR::performPCA(dataset)
  dataset = r_pca$dataset_ncp
} else {
  message("PCA is NOT suitable")
}

head(dataset)

# Second clustering and optimal k

r_atsc = evaluomeR::ATSC(data=dataset, alpha=0.1, k.range=c(3,10), cbi="kmeans", gold_standard=gold_standard_vector)

r_atsc$optimalK
r_atsc$trimmedRows
r_atsc$trimmedColumns
r_atsc$gold_standard_trimmed
r_atsc$trimmmedDataset

#evaluomeR::getRSKCAlpha(dataset, k=3, L1=3, max_alpha = 0.05)

