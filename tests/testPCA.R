library(evaluomeR)
library(RSKC)
library(sparcl)
library(parallel)

data("ontMetricsOBO")
data("golub")
data("nci60_k8")
data("ontMetrics")
data("bioMetrics")
data("rnaMetrics")
data("breastCancer")

seed=100

#Función auxiliar para series de ejecuciones
ejecutar_experimento_ATSC_sinPCA <- function(dataset, k.range, cbi="kmeans", seed=100, nEjec=5){
  tiempos <- numeric(nEjec)
  k_optimos <- numeric(nEjec)
  estabilidades <- numeric(nEjec)
  calidades <- numeric(nEjec)
  for(i in 1:nEjec){
    t <- system.time({
      r <- ATSC(dataset, k.range = k.range, cbi = cbi, seed = seed)
    })
    tiempos[i] <- t["elapsed"]
    k_optimos[i] <- r$optimalK_ATSC
    indice <- k_optimos[i] - min(k.range) + 1
    estabilidades[i] <- r$stab_ATSC[[indice]]
    calidades[i] <- r$qual_ATSC[[indice]]
  }
  return(list(tiempos = tiempos, media_tiempos = mean(tiempos), k_optimos = k_optimos,
              media_k = mean(k_optimos), estabilidades = estabilidades, media_stability = mean(estabilidades),
              calidades = calidades, media_quality = mean(calidades)))
}

#Función auxiliar para series de ejecuciones
ejecutar_experimento_ATSC_conPCA <- function(dataset, k.range, cbi="kmeans", correlation_threshold, numCores=1, seed=100, nEjec=5){
  tiempos <- numeric(nEjec)
  k_optimos <- numeric(nEjec)
  estabilidades <- numeric(nEjec)
  calidades <- numeric(nEjec)
  for(i in 1:nEjec){
    t <- system.time({
      r_cleanDataset = evaluomeR::cleanDataset(dataset, correlation_threshold)
      dataset_clean =  r_cleanDataset$dataset
      R =  r_cleanDataset$R
      
      pca_suitability = evaluomeR::PCASuitability(R, sig_level = 0.05)
      
      if (pca_suitability$pca_suitable) {
        message("PCA is suitable")
        r_pca = evaluomeR::performPCA(dataset_clean)
        dataset_postPCA = r_pca$dataset_ncp
      }
      r <- ATSC(dataset_postPCA, k.range = k.range, cbi = cbi, seed = seed, numCores = numCores)
    })
    tiempos[i] <- t["elapsed"]
    k_optimos[i] <- r$optimalK_ATSC
    indice <- k_optimos[i] - min(k.range) + 1
    estabilidades[i] <- r$stab_ATSC[[indice]]
    calidades[i] <- r$qual_ATSC[[indice]]
  }
  return(list(tiempos = tiempos, media_tiempos = mean(tiempos), k_optimos = k_optimos,
              media_k = mean(k_optimos), estabilidades = estabilidades, media_stability = mean(estabilidades),
              calidades = calidades, media_quality = mean(calidades)))
}


# PRUEBAS SUITABILITY PCA

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


# breastCancer

breastCancer_clean <- breastCancer[, -2]
t_PCA_breastCancer <- system.time({
  r_cleanDataset = evaluomeR::cleanDataset(breastCancer_clean, correlation_threshold = 0.85)
  dataset =  r_cleanDataset$dataset
  R =  r_cleanDataset$R
  
  pca_suitability = evaluomeR::PCASuitability(R, sig_level = 0.05)
  
  if (pca_suitability$pca_suitable) {
    message("PCA is suitable")
    r_pca = evaluomeR::performPCA(dataset)
    breastCancer_postPCA = r_pca$dataset_ncp
  }
})

cat("Métricas breastCancer antes de PCA:", ncol(breastCancer_clean)-1, "\n")
cat("Métricas breastCancer después de PCA:", ncol(breastCancer_postPCA)-1, "\n")
cat("Tiempo PCA breastCancer:", t_PCA_breastCancer["elapsed"], "seg\n")



# TESTS ATSC Y PCA

# ontMetricsOBO

ATSC1_sinPCA <- ejecutar_experimento_ATSC_sinPCA(ontMetricsOBO, k.range=c(2,6), cbi="kmeans", seed=100)

ATSC1_conPCA <- ejecutar_experimento_ATSC_conPCA(ontMetricsOBO, k.range=c(2,6), cbi="kmeans", correlation_threshold=0.95, seed=100)

cat("Tiempo ATSC ontMetricsOBO sin PCA:", ATSC1_sinPCA$media_tiempos, "\n")
cat("K óptimo ATSC ontMetricsOBO sin PCA:", ATSC1_sinPCA$media_k, "\n")
cat("Estabilidad ATSC ontMetricsOBO sin PCA:", ATSC1_sinPCA$media_stability, "\n")
cat("Calidad ATSC ontMetricsOBO sin PCA:", ATSC1_sinPCA$media_quality, "\n")

cat("Tiempo ATSC ontMetricsOBO con PCA:", ATSC1_conPCA$media_tiempos, "\n")
cat("K óptimo ATSC ontMetricsOBO con PCA:", ATSC1_conPCA$media_k, "\n")
cat("Estabilidad ATSC ontMetricsOBO con PCA:", ATSC1_conPCA$media_stability, "\n")
cat("Calidad ATSC ontMetricsOBO con PCA:", ATSC1_conPCA$media_quality, "\n")


# ontMetrics

ATSC2_sinPCA <- ejecutar_experimento_ATSC_sinPCA(ontMetrics_df, k.range=c(2,4), cbi="kmeans", seed=100)

ATSC2_conPCA <- ejecutar_experimento_ATSC_conPCA(ontMetrics_df, k.range=c(2,4), cbi="kmeans", correlation_threshold=0.90, seed=100)

cat("Tiempo ATSC ontMetrics sin PCA:", ATSC2_sinPCA$media_tiempos, "\n")
cat("K óptimo ATSC ontMetrics sin PCA:", ATSC2_sinPCA$media_k, "\n")
cat("Estabilidad ATSC ontMetrics sin PCA:", ATSC2_sinPCA$media_stability, "\n")
cat("Calidad ATSC ontMetrics sin PCA:", ATSC2_sinPCA$media_quality, "\n")

cat("Tiempo ATSC ontMetrics con PCA:", ATSC2_conPCA$media_tiempos, "\n")
cat("K óptimo ATSC ontMetrics con PCA:", ATSC2_conPCA$media_k, "\n")
cat("Estabilidad ATSC ontMetrics con PCA:", ATSC2_conPCA$media_stability, "\n")
cat("Calidad ATSC ontMetrics con PCA:", ATSC2_conPCA$media_quality, "\n")


# breastCancer

ATSC3_sinPCA <- ejecutar_experimento_ATSC_sinPCA(breastCancer_clean, k.range=c(2,6), cbi="kmeans", seed=100)

ATSC3_conPCA <- ejecutar_experimento_ATSC_conPCA(breastCancer_clean, k.range=c(2,6), cbi="kmeans", correlation_threshold=0.85, seed=100)

cat("Tiempo ATSC breastCancer sin PCA:", ATSC3_sinPCA$media_tiempos, "\n")
cat("K óptimo ATSC breastCancer sin PCA:", ATSC3_sinPCA$media_k, "\n")
cat("Estabilidad ATSC breastCancer sin PCA:", ATSC3_sinPCA$media_stability, "\n")
cat("Calidad ATSC breastCancer sin PCA:", ATSC3_sinPCA$media_quality, "\n")

cat("Tiempo ATSC breastCancer con PCA:", ATSC3_conPCA$media_tiempos, "\n")
cat("K óptimo ATSC breastCancer con PCA:", ATSC3_conPCA$media_k, "\n")
cat("Estabilidad ATSC breastCancer con PCA:", ATSC3_conPCA$media_stability, "\n")
cat("Calidad ATSC breastCancer con PCA:", ATSC3_conPCA$media_quality, "\n")

