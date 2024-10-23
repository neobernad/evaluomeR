library(evaluomeR)
library(RSKC)

data("ontMetricsOBO")
dataset = ontMetricsOBO

data("golub")
dataset = golub

# First, data cleaning

r_cleanDataset = evaluomeR::cleanDataset(dataset, correlation_threshold = 1)
dataset =  r_cleanDataset$dataset

pca_suitability = evaluomeR::PCASuitability(r_cleanDataset$R, sig_level = 0.05)
print(pca_suitability$pca_suitable)

if (pca_suitability$pca_suitable) {
  message("PCA is suitable")
  r_pca = evaluomeR::performPCA(dataset)
  dataset = r_pca$dataset_ncp
}

head(dataset)

# Second clustering and optimal k


r_atsc = evaluomeR::ATSC(data=dataset, k.range=c(3,10), cbi="kmeans")

r_atsc$optimalK
r_atsc$trimmedColumns
new_dataset = r_atsc$trimmmedDataset

