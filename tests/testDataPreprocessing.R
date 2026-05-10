library(evaluomeR)

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


