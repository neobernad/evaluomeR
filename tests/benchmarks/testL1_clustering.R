library(evaluomeR)
library(RSKC)
library(sparcl)
library(parallel)

data("ontMetrics")
data("nci60_k8")

seed=100

#Función auxiliar para series de ejecuciones
ejecutar_experimento_L1 <- function(dataset, k = NULL, clustering=clustering, seed = 100, nEjec = 1){
  tiempos <- numeric(nEjec)
  valores_L1 <- numeric(nEjec)
  for(i in 1:nEjec){
    t <- system.time({
      if(clustering == "kmeans"){
        L1 <- getRSKCL1Boundry(dataset, k=k, clustering=clustering, seed=seed)
      } else if(clustering == "hierarchical"){
        L1 <- getRSKCL1Boundry(dataset, clustering=clustering, seed=seed)
      }
    })
    tiempos[i] <- t["elapsed"]
    valores_L1[i] <- L1
  }
  return(list(tiempos = tiempos, media_tiempos = mean(tiempos), media_L1 = mean(valores_L1)))
}

#Función auxiliar para series de ejecuciones
ejecutar_experimento_L1_ATSC <- function(dataset, k.range=c(2, 4), cbi = "kmeans", clustering=clustering, L1=L1, seed = 100, nEjec = 1){
  tiempos <- numeric(nEjec)
  k_optimos <- numeric(nEjec)
  estabilidades <- numeric(nEjec)
  calidades <- numeric(nEjec)
  
  for(i in 1:nEjec){
    t <- system.time({
      r <- ATSC(dataset, k.range=k.range, cbi=cbi, L1=L1, numCores=1, seed=seed)
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


# NÚMERO DE MÉTRICAS (DATASETS)

# Fijamos k = 8, L1 = 3 (en el caso del cálculo de alpha), max_alpha = 0.1

alpha_metrica1 <- system.time({
  r_alpha_metrica1 <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=1, seed=100)
})

L1k_metrica1 <- system.time({
  r_L1k_metrica1 <- getRSKCL1Boundry(ontMetrics, k=8, clustering="kmeans", seed=100)
})

L1h_metrica1 <- system.time({
  r_L1h_metrica1 <- getRSKCL1Boundry(ontMetrics, clustering="hierarchical", seed=100)
})


alpha_metrica2 <- system.time({
  r_alpha_metrica2 <- getRSKCAlpha(nci60_k8, k=8, L1=3, max_alpha=0.10, numCores=1, seed=100)
})

L1k_metrica2 <- system.time({
  r_L1k_metrica2 <- getRSKCL1Boundry(nci60_k8, k=8, clustering="kmeans", seed=100)
})

L1h_metrica2 <- system.time({
  r_L1h_metrica2 <- getRSKCL1Boundry(nci60_k8, clustering="hierarchical", seed=100)
})

cat("Tiempo cálculo alpha bound para ontMetrics (19 métricas):", alpha_metrica1["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound (k-means clustering) para ontMetrics (19 métricas):", L1k_metrica1["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound (hierarchical clustering) para ontMetrics (19 métricas):", L1h_metrica1["elapsed"], "seg\n")

cat("Tiempo cálculo alpha bound para nci60_k8 (200 métricas):", alpha_metrica2["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound (k-means clustering) para nci60_k8 (200 métricas):", L1k_metrica2["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound (hierarchical clustering) para nci60_k8 (200 métricas):", L1h_metrica2["elapsed"], "seg\n")


# VALOR DE K ENTRE CLUSTERING METHODS 

# Fijamos dataset = ontMetrics

L1k_2k_ontMetrics <- ejecutar_experimento_L1(ontMetrics, k=2, clustering="kmeans", seed=100)

L1k_4k_ontMetrics <- ejecutar_experimento_L1(ontMetrics, k=4, clustering="kmeans", seed=100)

L1k_8k_ontMetrics <- ejecutar_experimento_L1(ontMetrics, k=8, clustering="kmeans", seed=100)

L1k_15k_ontMetrics <- ejecutar_experimento_L1(ontMetrics, k=15, clustering="kmeans", seed=100)

L1h_ontMetrics <- ejecutar_experimento_L1(ontMetrics, clustering="hierarchical", seed=100)


cat("Tiempo cálculo sparsity bound ontMetrics (k-means clustering, k=2):", L1k_2k_ontMetrics$media_tiempos, "seg\n")
cat("Tiempo cálculo sparsity bound ontMetrics (k-means clustering, k=4):", L1k_4k_ontMetrics$media_tiempos, "seg\n")
cat("Tiempo cálculo sparsity bound ontMetrics (k-means clustering, k=8):", L1k_8k_ontMetrics$media_tiempos, "seg\n")
cat("Tiempo cálculo sparsity bound ontMetrics (k-means clustering, k=15):", L1k_15k_ontMetrics$media_tiempos, "seg\n")
cat("Tiempo cálculo sparsity bound ontMetrics (hierarchical clustering):", L1h_ontMetrics$media_tiempos, "seg\n")

cat("Valor sparsity bound ontMetrics (k-means clustering, k=2):", L1k_2k_ontMetrics$media_L1, "\n")
cat("Valor sparsity bound ontMetrics (k-means clustering, k=4):", L1k_4k_ontMetrics$media_L1, "\n")
cat("Valor sparsity bound ontMetrics (k-means clustering, k=8):", L1k_8k_ontMetrics$media_L1, "\n")
cat("Valor sparsity bound ontMetrics (k-means clustering, k=15):", L1k_15k_ontMetrics$media_L1, "\n")
cat("Valor sparsity bound ontMetrics (hierarchical clustering):", L1h_ontMetrics$media_L1, "\n")


# Fijamos dataset = nci60_k8

L1k_2k_nci <- ejecutar_experimento_L1(nci60_k8, k=2, clustering="kmeans", seed=100)

L1k_4k_nci <- ejecutar_experimento_L1(nci60_k8, k=4, clustering="kmeans", seed=100)

L1k_8k_nci <- ejecutar_experimento_L1(nci60_k8, k=8, clustering="kmeans", seed=100)

L1k_15k_nci <- ejecutar_experimento_L1(nci60_k8, k=15, clustering="kmeans", seed=100)

L1h_nci <- ejecutar_experimento_L1(nci60_k8, clustering="hierarchical", seed=100)


cat("Tiempo cálculo sparsity bound nci60_k8 (k-means clustering, k=2):", L1k_2k_nci$media_tiempos, "seg\n")
cat("Tiempo cálculo sparsity bound nci60_k8 (k-means clustering, k=4):", L1k_4k_nci$media_tiempos, "seg\n")
cat("Tiempo cálculo sparsity bound nci60_k8 (k-means clustering, k=8):", L1k_8k_nci$media_tiempos, "seg\n")
cat("Tiempo cálculo sparsity bound nci60_k8 (k-means clustering, k=15):", L1k_15k_nci$media_tiempos, "seg\n")
cat("Tiempo cálculo sparsity bound nci60_k8 (hierarchical clustering):", L1h_nci$media_tiempos, "seg\n")

cat("Valor sparsity bound nci60_k8 (k-means clustering, k=2):", L1k_2k_nci$media_L1, "\n")
cat("Valor sparsity bound nci60_k8 (k-means clustering, k=4):", L1k_4k_nci$media_L1, "\n")
cat("Valor sparsity bound nci60_k8 (k-means clustering, k=8):", L1k_8k_nci$media_L1, "\n")
cat("Valor sparsity bound nci60_k8 (k-means clustering, k=15):", L1k_15k_nci$media_L1, "\n")
cat("Valor sparsity bound nci60_k8 (hierarchical clustering):", L1h_nci$media_L1, "\n")


# COMPARACIÓN VALOR K ÓPTIMO POST ATSC SEGÚN MÉTODO DE CLUSTERING DEL SPARSITY BOUND

clustering2k_ontMetrics <- ejecutar_experimento_L1_ATSC(ontMetrics, L1=2, seed=100)

clustering2k_nci60_k8 <- ejecutar_experimento_L1_ATSC(nci60_k8, L1=L1k_2k_nci$media_L1, seed=100)

clustering4k_nci60_k8 <- ejecutar_experimento_L1_ATSC(nci60_k8, L1=L1k_4k_nci$media_L1, seed=100)

clustering8k_nci60_k8 <- ejecutar_experimento_L1_ATSC(nci60_k8, L1=L1k_8k_nci$media_L1, seed=100)

clustering15k_nci60_k8 <- ejecutar_experimento_L1_ATSC(nci60_k8, L1=L1k_15k_nci$media_L1, seed=100)

clusteringh_nci60_k8 <- ejecutar_experimento_L1_ATSC(nci60_k8, L1=L1h_nci$media_L1, seed=100)

cat("Valor k-optimo post ATSC ontMetrics:", clustering2k_ontMetrics$media_k, "\n")
cat("Valor k-optimo post ATSC nci60_k8 (k-means sparsity clustering, k=2):", clustering2k_nci60_k8$media_k, "\n")
cat("Valor k-optimo post ATSC nci60_k8 (k-means sparsity clustering, k=4):", clustering4k_nci60_k8$media_k, "\n")
cat("Valor k-optimo post ATSC nci60_k8 (k-means sparsity clustering, k=8):", clustering8k_nci60_k8$media_k, "\n")
cat("Valor k-optimo post ATSC nci60_k8 (k-means sparsity clustering, k=15):", clustering15k_nci60_k8$media_k, "\n")
cat("Valor k-optimo post ATSC nci60_k8 (hierarchical sparsity clustering):", clusteringh_nci60_k8$media_k, "\n")

cat("Estabilidad post ATSC ontMetrics:", clustering2k_ontMetrics$media_stability, "\n")
cat("Estabilidad post ATSC nci60_k8 (k-means sparsity clustering, k=2):", clustering2k_nci60_k8$media_stability, "\n")
cat("Estabilidad post ATSC nci60_k8 (k-means sparsity clustering, k=4):", clustering4k_nci60_k8$media_stability, "\n")
cat("Estabilidad post ATSC nci60_k8 (k-means sparsity clustering, k=8):", clustering8k_nci60_k8$media_stability, "\n")
cat("Estabilidad post ATSC nci60_k8 (k-means sparsity clustering, k=15):", clustering15k_nci60_k8$media_stability, "\n")
cat("Estabilidad post ATSC nci60_k8 (hierarchical sparsity clustering):", clusteringh_nci60_k8$media_stability, "\n")

cat("Calidad post ATSC ontMetrics:", clustering2k_ontMetrics$media_quality, "\n")
cat("Calidad post ATSC nci60_k8 (k-means sparsity clustering, k=2):", clustering2k_nci60_k8$media_quality, "\n")
cat("Calidad post ATSC nci60_k8 (k-means sparsity clustering, k=4):", clustering4k_nci60_k8$media_quality, "\n")
cat("Calidad post ATSC nci60_k8 (k-means sparsity clustering, k=8):", clustering8k_nci60_k8$media_quality, "\n")
cat("Calidad post ATSC nci60_k8 (k-means sparsity clustering, k=15):", clustering15k_nci60_k8$media_quality, "\n")
cat("Calidad post ATSC nci60_k8 (hierarchical sparsity clustering):", clusteringh_nci60_k8$media_quality, "\n")

