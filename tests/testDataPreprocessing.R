library(evaluomeR)
library(RSKC)

data("ontMetricsOBO")
dataset = ontMetricsOBO

r_cleanDataset = evaluomeR::cleanDataset(dataset, correlation_threshold = 1)
dataset =  r_cleanDataset$dataset
R =  r_cleanDataset$R


pca_suitability = evaluomeR::PCASuitability(R, sig_level = 0.05)

if (pca_suitability$pca_suitable) {
  message("PCA is suitable")
  r_pca = evaluomeR::performPCA(dataset)
  #evaluomeR::plotPCA_fviz_screeplot(r_pca$pca)
  #evaluomeR::plotPCA_fviz_biplot(r_pca$pca)
  dataset = r_pca$dataset_ncp
}

head(dataset)



test = RSKC(dataset[-1], L=2, alpha=0.8, ncl=3)
union_vector = c(test$oE,test$oW)
print(union_vector)
length(unique(union_vector))

dataset_cleaned <- dataset[-union_vector, ]

