library(evaluomeR)
library(RSKC)

data("golub")
dataset = golub

# First, data cleaning

r_cleanDataset = evaluomeR::cleanDataset(dataset, correlation_threshold = 1)
dataset =  r_cleanDataset$dataset

pca_suitability = evaluomeR::PCASuitability(r_cleanDataset$R, sig_level = 0.05)
print(pca_suitability$pca_suitable)

#if (pca_suitability$pca_suitable) {
if (FALSE) {
  message("PCA is suitable")
  r_pca = evaluomeR::performPCA(dataset)
  dataset = r_pca$dataset_ncp
}

head(dataset)

# Second clustering and optimal k


r_atsc = evaluomeR::ATSC(data=dataset, k.range=c(3,10), cbi="clara")

r_atsc$trimmedRows
r_atsc$trimmedColumns
r_atsc = r_atsc$trimmmedDataset

