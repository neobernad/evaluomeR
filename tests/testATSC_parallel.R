library(evaluomeR)
library(RSKC)
library(sparcl)
library(parallel)

data("ontMetrics")
data("golub")
data("nci60_k8")
data("breastCancer")

seed=100

#Función auxiliar para series de ejecuciones
ATSC_partes <- function(data, k.range=c(2,15), bs=100, cbi="kmeans", max_alpha = 0.1,
                                L1=NULL, alpha=NULL, gold_standard=NULL,
                                seed=NULL, numCores=1, clusteringSparsity="kmeans", nEjec=5) {
  tiempos_df <- data.frame(
    t1_stability = numeric(nEjec),
    t1_quality = numeric(nEjec),
    t2_optimalK = numeric(nEjec),
    t3_L1 = numeric(nEjec),
    t4_alpha = numeric(nEjec),
    t5_rskc = numeric(nEjec),
    t6_stability = numeric(nEjec),
    t6_quality = numeric(nEjec),
    t7_optimalK = numeric(nEjec),
    total = numeric(nEjec)
  )
  
  for(i in 1:nEjec){
    L1 = NULL
    alpha=NULL
    all_metrics = TRUE
    
    message(paste0("Computing optimal k value with '", cbi, "'"))
    data = as.data.frame(SummarizedExperiment::assay(data))
    t1_stability <- system.time({
      stabRange = evaluomeR::stabilityRange(data=data, cbi=cbi, k=k.range, bs=bs,
                                            all_metrics = TRUE,
                                            gold_standard=gold_standard, seed=seed, numCores=1)
      stab = evaluomeR::standardizeStabilityData(stabRange, k.range = k.range)
    })
    
    t1_quality <- system.time({
      qualRange = evaluomeR::qualityRange(data, k=k.range, cbi=cbi, all_metrics = TRUE, seed=seed, numCores=1)
      qual = evaluomeR::standardizeQualityData(qualRange, k.range = k.range)
    })
    
    tiempos_df$t1_stability[i] <- t1_stability["elapsed"]
    tiempos_df$t1_quality[i] <- t1_quality["elapsed"]
    
    t2_optimalK <- system.time({
      rOptimalK = evaluomeR::getOptimalKValue(stabRange, qualRange)
      optimalK = as.numeric(rOptimalK$Global_optimal_k)
      message(paste0("Optimal k: ", optimalK))
    })
    
    tiempos_df$t2_optimalK[i] <- t2_optimalK["elapsed"]
    
    message("Determining best L1 and alpha parameter automatically, it might take a while...")
    
    t3_L1 <- system.time({
      L1 = evaluomeR::getRSKCL1Boundry(data, k=optimalK, clustering=clusteringSparsity, seed=seed)
    })
    
    tiempos_df$t3_L1[i] <- t3_L1["elapsed"]
    
    t4_alpha <- system.time({
      best_alphas = evaluomeR::getRSKCAlpha(data, k=optimalK, L1=L1,
                                            max_alpha = max_alpha, seed=seed, numCores=numCores)
      alpha = best_alphas$best_alpha_qual
    })
    
    tiempos_df$t4_alpha[i] <- t4_alpha["elapsed"]
    
    message(paste0("\tUsing L1 '", L1, "' and alpha '", alpha, "'"))
    
    t5_rskc <- system.time({
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
    })
    
    tiempos_df$t5_rskc[i] <- t5_rskc["elapsed"]
    t6_stability <- system.time({
      stabRange_ATSC = evaluomeR::stabilityRange(data=data_trimmed, cbi=cbi, k=k.range, bs=bs,
                                                 all_metrics = TRUE,
                                                 gold_standard=gold_standard, seed=seed, numCores=1)
      stab_ATSC = evaluomeR::standardizeStabilityData(stabRange_ATSC, k.range = k.range)
    })
    
    t6_quality <- system.time({
      qualRange_ATSC = evaluomeR::qualityRange(data=data_trimmed, cbi=cbi, k=k.range,
                                               all_metrics = TRUE, seed=seed, numCores=1)
      qual_ATSC = evaluomeR::standardizeQualityData(qualRange_ATSC, k.range = k.range)
    })

    tiempos_df$t6_stability[i] <- t6_stability["elapsed"]
    tiempos_df$t6_quality[i] <- t6_quality["elapsed"] 
    
    t7_optimalK <- system.time({
      rOptimalK_ATSC = evaluomeR::getOptimalKValue(stabRange_ATSC, qualRange_ATSC)
      optimalK_ATSC = as.numeric(rOptimalK_ATSC$Global_optimal_k)
      message(paste0("New optimal k with ATSC: ", optimalK_ATSC))
    })
    
    tiempos_df$t7_optimalK[i] <- t7_optimalK["elapsed"]
    tiempos_df$total[i] <- sum(tiempos_df[i, 1:9])
  }
  medias <- as.list(colMeans(tiempos_df))
  return(list(tiempos = tiempos_df, media = medias))
  
}

# ontMetrics

ATSC_ontMetrics_sec <- ATSC_partes(ontMetrics, k.range=c(3,10), cbi="kmeans", max_alpha=0.10, numCores=1, seed=100)

cat("Fase 1: Indice de estabilidad", ATSC_ontMetrics_sec$media$t1_stability, "seg\n")
cat("Fase 1: Indice de calidad", ATSC_ontMetrics_sec$media$t1_quality, "seg\n")
cat("Fase 2: K óptimo", ATSC_ontMetrics_sec$media$t2_optimalK, "seg\n")
cat("Fase 3: L1", ATSC_ontMetrics_sec$media$t3_L1, "seg\n")
cat("Fase 4: Alpha óptimo", ATSC_ontMetrics_sec$media$t4_alpha, "seg\n")
cat("Fase 5: RSKC", ATSC_ontMetrics_sec$media$t5_rskc, "seg\n")
cat("Fase 6: Indice de estabilidad", ATSC_ontMetrics_sec$media$t6_stability, "seg\n")
cat("Fase 6: Indice de calidad", ATSC_ontMetrics_sec$media$t6_quality, "seg\n")
cat("Tiempo ATSC total", ATSC_ontMetrics_sec$media$total, "seg\n")

ATSC_ontMetrics_par <- ATSC_partes(ontMetrics, k.range=c(3,10), cbi="kmeans", max_alpha=0.10, numCores=12, seed=100)

cat("Fase 1: Indice de estabilidad", ATSC_ontMetrics_par$media$t1_stability, "seg\n")
cat("Fase 1: Indice de calidad", ATSC_ontMetrics_par$media$t1_quality, "seg\n")
cat("Fase 2: K óptimo", ATSC_ontMetrics_par$media$t2_optimalK, "seg\n")
cat("Fase 3: L1", ATSC_ontMetrics_par$media$t3_L1, "seg\n")
cat("Fase 4: Alpha óptimo", ATSC_ontMetrics_par$media$t4_alpha, "seg\n")
cat("Fase 5: RSKC", ATSC_ontMetrics_par$media$t5_rskc, "seg\n")
cat("Fase 6: Indice de estabilidad", ATSC_ontMetrics_par$media$t6_stability, "seg\n")
cat("Fase 6: Indice de calidad", ATSC_ontMetrics_par$media$t6_quality, "seg\n")
cat("Tiempo ATSC paralelo total", ATSC_ontMetrics_par$media$total, "seg\n")


# golub
r_cleanDataset = cleanDataset(golub, correlation_threshold = 1)
golub =  r_cleanDataset$dataset
ATSC_golub_sec <- ATSC_partes(golub, k.range=c(3,10), cbi="kmeans", max_alpha=0.10, numCores=1, seed=100)

cat("Fase 1: Indice de estabilidad", ATSC_golub_sec$media$t1_stability, "seg\n")
cat("Fase 1: Indice de calidad", ATSC_golub_sec$media$t1_quality, "seg\n")
cat("Fase 2: K óptimo", ATSC_golub_sec$media$t2_optimalK, "seg\n")
cat("Fase 3: L1", ATSC_golub_sec$media$t3_L1, "seg\n")
cat("Fase 4: Alpha óptimo", ATSC_golub_sec$media$t4_alpha, "seg\n")
cat("Fase 5: RSKC", ATSC_golub_sec$media$t5_rskc, "seg\n")
cat("Fase 6: Indice de estabilidad", ATSC_golub_sec$media$t6_stability, "seg\n")
cat("Fase 6: Indice de calidad", ATSC_golub_sec$media$t6_quality, "seg\n")
cat("Tiempo ATSC total", ATSC_golub_sec$media$total, "seg\n")

ATSC_golub_par <- ATSC_partes(golub, k.range=c(3,10), cbi="kmeans", max_alpha=0.10, numCores=12, seed=100)

cat("Fase 1: Indice de estabilidad", ATSC_golub_par$media$t1_stability, "seg\n")
cat("Fase 1: Indice de calidad", ATSC_golub_par$media$t1_quality, "seg\n")
cat("Fase 2: K óptimo", ATSC_golub_par$media$t2_optimalK, "seg\n")
cat("Fase 3: L1", ATSC_golub_par$media$t3_L1, "seg\n")
cat("Fase 4: Alpha óptimo", ATSC_golub_par$media$t4_alpha, "seg\n")
cat("Fase 5: RSKC", ATSC_golub_par$media$t5_rskc, "seg\n")
cat("Fase 6: Indice de estabilidad", ATSC_golub_par$media$t6_stability, "seg\n")
cat("Fase 6: Indice de calidad", ATSC_golub_par$media$t6_quality, "seg\n")
cat("Tiempo ATSC paralelo total", ATSC_golub_par$media$total, "seg\n")


# nci60_k8

ATSC_nci60_k8_sec <- ATSC_partes(nci60_k8, k.range=c(3,10), cbi="kmeans", max_alpha=0.10, numCores=1, seed=100)

cat("Fase 1: Indice de estabilidad", ATSC_nci60_k8_sec$media$t1_stability, "seg\n")
cat("Fase 1: Indice de calidad", ATSC_nci60_k8_sec$media$t1_quality, "seg\n")
cat("Fase 2: K óptimo", ATSC_nci60_k8_sec$media$t2_optimalK, "seg\n")
cat("Fase 3: L1", ATSC_nci60_k8_sec$media$t3_L1, "seg\n")
cat("Fase 4: Alpha óptimo", ATSC_nci60_k8_sec$media$t4_alpha, "seg\n")
cat("Fase 5: RSKC", ATSC_nci60_k8_sec$media$t5_rskc, "seg\n")
cat("Fase 6: Indice de estabilidad", ATSC_nci60_k8_sec$media$t6_stability, "seg\n")
cat("Fase 6: Indice de calidad", ATSC_nci60_k8_sec$media$t6_quality, "seg\n")
cat("Tiempo ATSC total", ATSC_nci60_k8_sec$media$total, "seg\n")

ATSC_nci60_k8_par <- ATSC_partes(nci60_k8, k.range=c(3,10), cbi="kmeans", max_alpha=0.10, numCores=12, seed=100)

cat("Fase 1: Indice de estabilidad", ATSC_nci60_k8_par$media$t1_stability, "seg\n")
cat("Fase 1: Indice de calidad", ATSC_nci60_k8_par$media$t1_quality, "seg\n")
cat("Fase 2: K óptimo", ATSC_nci60_k8_par$media$t2_optimalK, "seg\n")
cat("Fase 3: L1", ATSC_nci60_k8_par$media$t3_L1, "seg\n")
cat("Fase 4: Alpha óptimo", ATSC_nci60_k8_par$media$t4_alpha, "seg\n")
cat("Fase 5: RSKC", ATSC_nci60_k8_par$media$t5_rskc, "seg\n")
cat("Fase 6: Indice de estabilidad", ATSC_nci60_k8_par$media$t6_stability, "seg\n")
cat("Fase 6: Indice de calidad", ATSC_nci60_k8_par$media$t6_quality, "seg\n")
cat("Tiempo ATSC paralelo total", ATSC_nci60_k8_par$media$total, "seg\n")


# breastCancer

ATSC_breastCancer_sec <- ATSC_partes(breastCancer, k.range=c(3,10), cbi="kmeans", max_alpha=0.10, numCores=1, seed=100)

cat("Fase 1: Indice de estabilidad", ATSC_breastCancer_sec$media$t1_stability, "seg\n")
cat("Fase 1: Indice de calidad", ATSC_breastCancer_sec$media$t1_quality, "seg\n")
cat("Fase 2: K óptimo", ATSC_breastCancer_sec$media$t2_optimalK, "seg\n")
cat("Fase 3: L1", ATSC_breastCancer_sec$media$t3_L1, "seg\n")
cat("Fase 4: Alpha óptimo", ATSC_breastCancer_sec$media$t4_alpha, "seg\n")
cat("Fase 5: RSKC", ATSC_breastCancer_sec$media$t5_rskc, "seg\n")
cat("Fase 6: Indice de estabilidad", ATSC_breastCancer_sec$media$t6_stability, "seg\n")
cat("Fase 6: Indice de calidad", ATSC_breastCancer_sec$media$t6_quality, "seg\n")
cat("Tiempo ATSC total", ATSC_breastCancer_sec$media$total, "seg\n")

ATSC_breastCancer_par <- ATSC_partes(breastCancer, k.range=c(3,10), cbi="kmeans", max_alpha=0.10, numCores=12, seed=100)

cat("Fase 1: Indice de estabilidad", ATSC_breastCancer_par$media$t1_stability, "seg\n")
cat("Fase 1: Indice de calidad", ATSC_breastCancer_par$media$t1_quality, "seg\n")
cat("Fase 2: K óptimo", ATSC_breastCancer_par$media$t2_optimalK, "seg\n")
cat("Fase 3: L1", ATSC_breastCancer_par$media$t3_L1, "seg\n")
cat("Fase 4: Alpha óptimo", ATSC_breastCancer_par$media$t4_alpha, "seg\n")
cat("Fase 5: RSKC", ATSC_breastCancer_par$media$t5_rskc, "seg\n")
cat("Fase 6: Indice de estabilidad", ATSC_breastCancer_par$media$t6_stability, "seg\n")
cat("Fase 6: Indice de calidad", ATSC_breastCancer_par$media$t6_quality, "seg\n")
cat("Tiempo ATSC paralelo total", ATSC_breastCancer_par$media$total, "seg\n")


ATSC_breastCancer_par_PCA <- ejecutar_experimento_ATSC_conPCA(breastCancer, k.range=c(2,6), cbi="kmeans", correlation_threshold=0.85, seed=100, numCores=12)

cat("Tiempo ATSC paralelo breastCancer con PCA:", ATSC_breastCancer_par_PCA$media_tiempos, "seg\n")

