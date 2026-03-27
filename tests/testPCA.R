library(evaluomeR)

data("ontMetricsOBO")
data("golub")
data("nci60_k8")
data("ontMetrics")
data("bioMetrics")
data("rnaMetrics")


# ontMetricsOBO

t_PCA_OBO <- system.time({
  r_cleanDataset = evaluomeR::cleanDataset(ontMetricsOBO, correlation_threshold = 0.95)
  dataset =  r_cleanDataset$dataset
  R =  r_cleanDataset$R
  
  pca_suitability = evaluomeR::PCASuitability(R, sig_level = 0.05)
  
  if (pca_suitability$pca_suitable) {
    message("PCA is suitable")
    r_pca = evaluomeR::performPCA(dataset)
    ontMetricsOBO_postPCA = r_pca$dataset_ncp
  }
})

cat("Métricas ontMetricsOBO antes de PCA:", ncol(ontMetricsOBO)-1, "\n")
cat("Métricas ontMetricsOBO después de PCA:", ncol(ontMetricsOBO_postPCA)-1, "\n")
cat("Tiempo PCA ontMetricsOBO:", t_PCA_OBO["elapsed"], "seg\n")


# ontMetrics

if (is(ontMetrics, "SummarizedExperiment")) {
  df <- as.data.frame(SummarizedExperiment::assay(ontMetrics))
  ontMetrics_df <- df
}

t_PCA_ontMetrics <- system.time({
  r_cleanDataset = evaluomeR::cleanDataset(ontMetrics_df, correlation_threshold = 0.95)
  dataset =  r_cleanDataset$dataset
  R =  r_cleanDataset$R
  
  pca_suitability = evaluomeR::PCASuitability(R, sig_level = 0.05)
  
  if (pca_suitability$pca_suitable) {
    message("PCA is suitable")
    r_pca = evaluomeR::performPCA(dataset)
    ontMetrics_postPCA = r_pca$dataset_ncp
  }
})

cat("Métricas ontMetrics antes de PCA:", ncol(ontMetrics_df)-1, "\n")
cat("Métricas ontMetrics después de PCA:", ncol(ontMetrics_postPCA)-1, "\n")
cat("Tiempo PCA ontMetrics:", t_PCA_ontMetrics["elapsed"], "seg\n")


# nci60_k8

t_PCA_nci60_k8 <- system.time({
  r_cleanDataset = evaluomeR::cleanDataset(nci60_k8, correlation_threshold = 0.50)
  dataset =  r_cleanDataset$dataset
  R =  r_cleanDataset$R
  
  pca_suitability = evaluomeR::PCASuitability(R, sig_level = 0.05)
  
  if (pca_suitability$pca_suitable) {
    message("PCA is suitable")
    r_pca = evaluomeR::performPCA(dataset)
    nci60_k8_postPCA = r_pca$dataset_ncp
  }
})

cat("Métricas nci60_k8 antes de PCA:", ncol(nci60_k8)-1, "\n")
cat("Métricas nci60_k8 después de PCA:", ncol(nci60_k8_postPCA)-1, "\n")
cat("Tiempo PCA nci60_k8:", t_PCA_nci60_k8["elapsed"], "seg\n")


# nci60_k10

t_PCA_nci60_k10 <- system.time({
  r_cleanDataset = evaluomeR::cleanDataset(nci60_k10, correlation_threshold = 0.50)
  dataset =  r_cleanDataset$dataset
  R =  r_cleanDataset$R
  
  pca_suitability = evaluomeR::PCASuitability(R, sig_level = 0.05)
  
  if (pca_suitability$pca_suitable) {
    message("PCA is suitable")
    r_pca = evaluomeR::performPCA(dataset)
    nci60_k10_postPCA = r_pca$dataset_ncp
  }
})

cat("Métricas nci60_k10 antes de PCA:", ncol(nci60_k10)-1, "\n")
cat("Métricas nci60_k10 después de PCA:", ncol(nci60_k10_postPCA)-1, "\n")
cat("Tiempo PCA nci60_k10:", t_PCA_nci60_k10["elapsed"], "seg\n")


# golub

t_PCA_golub <- system.time({
  r_cleanDataset = evaluomeR::cleanDataset(golub, correlation_threshold = 0.65)
  dataset =  r_cleanDataset$dataset
  R =  r_cleanDataset$R
  
  pca_suitability = evaluomeR::PCASuitability(R, sig_level = 0.05)
  
  if (pca_suitability$pca_suitable) {
    message("PCA is suitable")
    r_pca = evaluomeR::performPCA(dataset)
    golub_postPCA = r_pca$dataset_ncp
  }
})

cat("Métricas golub antes de PCA:", ncol(golub)-1, "\n")
cat("Métricas golub después de PCA:", ncol(golub_postPCA)-1, "\n")
cat("Tiempo PCA golub:", t_PCA_golub["elapsed"], "seg\n")


# TESTS ATSC Y PCA

# ontMetricsOBO

t_ATSC1_sinPCA <- system.time({
  r_ontMetricsOBO = ATSC(ontMetricsOBO, k.range=c(2,6), cbi="kmeans")
})

t_ATSC1_conPCA <- system.time({
  r_cleanDataset = evaluomeR::cleanDataset(ontMetricsOBO, correlation_threshold = 0.95)
  dataset =  r_cleanDataset$dataset
  R =  r_cleanDataset$R
  
  pca_suitability = evaluomeR::PCASuitability(R, sig_level = 0.05)
  
  if (pca_suitability$pca_suitable) {
    message("PCA is suitable")
    r_pca = evaluomeR::performPCA(dataset)
    ontMetricsOBO_postPCA = r_pca$dataset_ncp
  }
  r_ontMetricsOBO_postPCA = ATSC(ontMetricsOBO_postPCA, k.range=c(2,6), cbi="kmeans")
})

cat("Tiempo ATSC ontMetricsOBO sin PCA:", t_ATSC1_sinPCA["elapsed"], "seg\n")
cat("Tiempo ATSC ontMetricsOBO con PCA:", t_ATSC1_conPCA["elapsed"], "seg\n")


# ontMetrics

t_ATSC2_sinPCA <- system.time({
  r_rnaMetrics = ATSC(ontMetrics, k.range=c(2,4), cbi="kmeans")
})

t_ATSC2_conPCA <- system.time({
  r_cleanDataset = evaluomeR::cleanDataset(ontMetrics_df, correlation_threshold = 0.90)
  dataset =  r_cleanDataset$dataset
  R =  r_cleanDataset$R
  
  pca_suitability = evaluomeR::PCASuitability(R, sig_level = 0.05)
  
  if (pca_suitability$pca_suitable) {
    message("PCA is suitable")
    r_pca = evaluomeR::performPCA(dataset)
    ontMetrics_postPCA = r_pca$dataset_ncp
  }
  r_ontMetrics_postPCA = ATSC(ontMetrics_postPCA, k.range=c(2,4), cbi="kmeans")
})

cat("Tiempo ATSC ontMetrics sin PCA:", t_ATSC2_sinPCA["elapsed"], "seg\n")
cat("Tiempo ATSC ontMetrics con PCA:", t_ATSC2_conPCA["elapsed"], "seg\n")
