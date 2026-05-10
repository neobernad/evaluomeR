library(evaluomeR)
library(RSKC)
library(sparcl)
library(parallel)

data("ontMetrics")
data("rnaMetrics")
data("golub")
data("nci60_k8")

seed=100

#Función auxiliar para series de ejecuciones
ejecutar_experimento_Alpha <- function(dataset, k, L1, max_alpha, numCores, seed, nEjec){
  tiempos <- numeric(nEjec)
  for(i in 1:nEjec){
    t <- system.time({
      r <- getRSKCAlpha(dataset, k=k, L1=L1, max_alpha=max_alpha, numCores=numCores, seed=seed)
    })
    tiempos[i] <- t["elapsed"]
  }
  return(list(tiempos = tiempos, media = mean(tiempos)))
}


# Tests para comprobar parámetro a parámetro la diferencia de rendimiento entre versiones


# NÚMERO DE MÉTRICAS (DATASETS)

# Fijamos k = 8, L1 = 3, max_alpha = 0.1 y numCores = 11

alpha_metrica1_sec <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=1, seed=100, nEjec=5)

alpha_metrica1_par <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=11, seed=100, nEjec=5)


alpha_metrica2_sec <- ejecutar_experimento_Alpha(nci60_k8, k=8, L1=3, max_alpha=0.10, numCores=1, seed=100, nEjec=5)

alpha_metrica2_par <- ejecutar_experimento_Alpha(nci60_k8, k=8, L1=3, max_alpha=0.10, numCores=11, seed=100, nEjec=5)

cat("Tiempo ontMetrics sec", alpha_metrica1_sec$media, "seg\n")
cat("Tiempo nci60_k8 sec", alpha_metrica2_sec$media, "seg\n")

cat("Tiempo ontMetrics par", alpha_metrica1_par$media, "seg\n")
cat("Tiempo nci60_k8 par", alpha_metrica2_par$media, "seg\n")


# NÚMERO DE CORES

# Fijamos dataset = ontMetrics, k = 8, L1 = 3 y max_alpha = 0.1

alpha_1cores <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=1, seed=100, nEjec=5)

alpha_2cores <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=2, seed=100, nEjec=5)

alpha_4cores <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=4, seed=100, nEjec=5)

alpha_6cores <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=6, seed=100, nEjec=5)

alpha_8cores <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=8, seed=100, nEjec=5)

alpha_11cores <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=11, seed=100, nEjec=5)

cat("Tiempo secuencial (1 core):", alpha_1cores$media, "seg\n")
cat("Tiempo 2 cores:", alpha_2cores$media, "seg\n")
cat("Tiempo 4 cores:", alpha_4cores$media, "seg\n")
cat("Tiempo 6 cores:", alpha_6cores$media, "seg\n")
cat("Tiempo 8 cores:", alpha_8cores$media, "seg\n")
cat("Tiempo 11 cores:", alpha_11cores$media, "seg\n")


# NÚMERO DE CLUSTERS (k)

# Fijamos dataset = ontMetrics, L1 = 3, max_alpha = 0.1, numCores = 11

alpha_2k_sec <- ejecutar_experimento_Alpha(ontMetrics, k=2, L1=3, max_alpha=0.10, numCores=1, seed=100, nEjec=5)

alpha_2k_par <- ejecutar_experimento_Alpha(ontMetrics, k=2, L1=3, max_alpha=0.10, numCores=11, seed=100, nEjec=5)

alpha_4k_sec <- ejecutar_experimento_Alpha(ontMetrics, k=4, L1=3, max_alpha=0.10, numCores=1, seed=100, nEjec=5)

alpha_4k_par <- ejecutar_experimento_Alpha(ontMetrics, k=4, L1=3, max_alpha=0.10, numCores=11, seed=100, nEjec=5)

alpha_8k_sec <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=1, seed=100, nEjec=5)

alpha_8k_par <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=11, seed=100, nEjec=5)

alpha_15k_sec <- ejecutar_experimento_Alpha(ontMetrics, k=15, L1=3, max_alpha=0.10, numCores=1, seed=100, nEjec=5)

alpha_15k_par <- ejecutar_experimento_Alpha(ontMetrics, k=15, L1=3, max_alpha=0.10, numCores=11, seed=100, nEjec=5)

cat("Tiempo k = 2 sec:", alpha_2k_sec$media, "seg\n")
cat("Tiempo k = 4 sec:", alpha_4k_sec$media, "seg\n")
cat("Tiempo k = 8 sec:", alpha_8k_sec$media, "seg\n")
cat("Tiempo k = 15 sec:", alpha_15k_sec$media, "seg\n")

cat("Tiempo k = 2 par:", alpha_2k_par$media, "seg\n")
cat("Tiempo k = 4 par:", alpha_4k_par$media, "seg\n")
cat("Tiempo k = 8 par:", alpha_8k_par$media, "seg\n")
cat("Tiempo k = 15 par:", alpha_15k_par$media, "seg\n")


# ALPHA MÁXIMO

# Fijamos dataset = ontMetrics, L1 = 3, k = 8 y numCores = 11

alpha_0.05_sec <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.05, numCores=1, seed=100, nEjec=5)

alpha_0.05_par <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.05, numCores=11, seed=100, nEjec=5)

alpha_0.07_sec <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.07, numCores=1, seed=100, nEjec=5)

alpha_0.07_par <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.07, numCores=11, seed=100, nEjec=5)

alpha_0.1_sec <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.1, numCores=1, seed=100, nEjec=5)

alpha_0.1_par <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.1, numCores=11, seed=100, nEjec=5)

alpha_0.15_sec <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.15, numCores=1, seed=100, nEjec=5)

alpha_0.15_par <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.15, numCores=11, seed=100, nEjec=5)

alpha_0.20_sec <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.2, numCores=1, seed=100, nEjec=5)

alpha_0.20_par <- ejecutar_experimento_Alpha(ontMetrics, k=8, L1=3, max_alpha=0.2, numCores=11, seed=100, nEjec=5)

cat("Tiempo max_alpha = 0.05 sec:", alpha_0.05_sec$media, "seg\n")
cat("Tiempo max_alpha = 0.07 sec:", alpha_0.07_sec$media, "seg\n")
cat("Tiempo max_alpha = 0.1 sec:", alpha_0.1_sec$media, "seg\n")
cat("Tiempo max_alpha = 0.15 sec:", alpha_0.15_sec$media, "seg\n")
cat("Tiempo max_alpha = 0.20 sec:", alpha_0.20_sec$media, "seg\n")

cat("Tiempo max_alpha = 0.05 par:", alpha_0.05_par$media, "seg\n")
cat("Tiempo max_alpha = 0.07 par:", alpha_0.07_par$media, "seg\n")
cat("Tiempo max_alpha = 0.1 par:", alpha_0.1_par$media, "seg\n")
cat("Tiempo max_alpha = 0.15 par:", alpha_0.15_par$media, "seg\n")
cat("Tiempo max_alpha = 0.20 par:", alpha_0.20_par$media, "seg\n")

