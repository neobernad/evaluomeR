library(evaluomeR)

data("ontMetricsOBO")
dataset = ontMetricsOBO

r_cleanDataset = evaluomeR::cleanDataset(dataset, correlation_threshold = 0.98)
dataset =  r_cleanDataset$dataset
R =  r_cleanDataset$R


pca_suitable = evaluomeR::PCASuitability(R, sig_level = 0.05)


nFactors = evaluomeR::determineNumberOfFactors(dataset)


r_pca = evaluomeR::performPCA(dataset = dataset, ncp = nFactors, scale = TRUE)
evaluomeR::plotPCA_fviz_screeplot(r_pca$pca)
evaluomeR::plotPCA_fviz_biplot(r_pca$pca)

head(r_pca$dataset_ncp)
