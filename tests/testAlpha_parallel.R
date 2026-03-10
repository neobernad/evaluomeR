library(evaluomeR)
library(RSKC)
library(sparcl)
library(parallel)

data("ontMetrics")
data("rnaMetrics")
data("golub")
data("nci60_k8")

seed=100

# Tests para comprobar parámetro a parámetro la diferencia de rendimiento entre versiones


# NÚMERO DE MÉTRICAS (DATASETS)

# Fijamos k = 8, L1 = 3, max_alpha = 0.1 y numCores = 6

alpha_metrica1_sec <- system.time({
  r_alpha_metrica1_sec <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=1)
})

alpha_metrica1_par <- system.time({
  r_alpha_metrica1_par <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=6)
})

alpha_metrica2_sec <- system.time({
  r_alpha_metrica2_sec <- getRSKCAlpha(nci60_k8, k=8, L1=3, max_alpha=0.10, numCores=1)
})

alpha_metrica2_par <- system.time({
  r_alpha_metrica2_par <- getRSKCAlpha(nci60_k8, k=8, L1=3, max_alpha=0.10, numCores=6)
})

cat("Tiempo ontMetrics sec", alpha_metrica1_sec["elapsed"], "seg\n")
cat("Tiempo nci60_k8 sec", alpha_metrica2_sec["elapsed"], "seg\n")

cat("Tiempo ontMetrics par", alpha_metrica1_par["elapsed"], "seg\n")
cat("Tiempo nci60_k8 par", alpha_metrica2_par["elapsed"], "seg\n")


# NÚMERO DE CORES

# Fijamos dataset = ontMetrics, k = 8, L1 = 3 y max_alpha = 0.1

alpha_1cores <- system.time({
  r_alpha_1cores <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=1)
})

alpha_2cores <- system.time({
  r_alpha_2cores <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=2)
})

alpha_4cores <- system.time({
  r_alpha_4cores <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=4)
})

alpha_6cores <- system.time({
  r_alpha_6cores <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=6)
})

alpha_8cores <- system.time({
  r_alpha_8cores <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=8)
})

alpha_11cores <- system.time({
  r_alpha_11cores <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=11)
})

cat("Tiempo secuencial (1 core):", alpha_1cores["elapsed"], "seg\n")
cat("Tiempo 2 cores:", alpha_2cores["elapsed"], "seg\n")
cat("Tiempo 4 cores:", alpha_4cores["elapsed"], "seg\n")
cat("Tiempo 6 cores:", alpha_6cores["elapsed"], "seg\n")
cat("Tiempo 8 cores:", alpha_8cores["elapsed"], "seg\n")
cat("Tiempo 11 cores:", alpha_11cores["elapsed"], "seg\n")


# NÚMERO DE CLUSTERS (k)

# Fijamos dataset = ontMetrics, L1 = 3, max_alpha = 0.1, numCores = 6

alpha_2k_sec <- system.time({
  r_alpha_2k_sec <- getRSKCAlpha(ontMetrics, k=2, L1=3, max_alpha=0.10, numCores=1)
})

alpha_2k_par <- system.time({
  r_alpha_2k_par <- getRSKCAlpha(ontMetrics, k=2, L1=3, max_alpha=0.10, numCores=6)
})

alpha_4k_sec <- system.time({
  r_alpha_4k_sec <- getRSKCAlpha(ontMetrics, k=4, L1=3, max_alpha=0.10, numCores=1)
})

alpha_4k_par <- system.time({
  r_alpha_4k_par <- getRSKCAlpha(ontMetrics, k=4, L1=3, max_alpha=0.10, numCores=6)
})

alpha_8k_sec <- system.time({
  r_alpha_8k_sec <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=1)
})

alpha_8k_par <- system.time({
  r_alpha_8k_par <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.10, numCores=6)
})

alpha_15k_sec <- system.time({
  r_alpha_15k_sec <- getRSKCAlpha(ontMetrics, k=15, L1=3, max_alpha=0.10, numCores=1)
})

alpha_15k_par <- system.time({
  r_alpha_15k_par <- getRSKCAlpha(ontMetrics, k=15, L1=3, max_alpha=0.10, numCores=6)
})

cat("Tiempo k = 2 sec:", alpha_2k_sec["elapsed"], "seg\n")
cat("Tiempo k = 4 sec:", alpha_4k_sec["elapsed"], "seg\n")
cat("Tiempo k = 8 sec:", alpha_8k_sec["elapsed"], "seg\n")
cat("Tiempo k = 15 sec:", alpha_15k_sec["elapsed"], "seg\n")

cat("Tiempo k = 2 par:", alpha_2k_par["elapsed"], "seg\n")
cat("Tiempo k = 4 par:", alpha_4k_par["elapsed"], "seg\n")
cat("Tiempo k = 8 par:", alpha_8k_par["elapsed"], "seg\n")
cat("Tiempo k = 15 par:", alpha_15k_par["elapsed"], "seg\n")


# ALPHA MÁXIMO

# Fijamos dataset = ontMetrics, L1 = 3, k = 8 y numCores = 6

alpha_0.05_sec <- system.time({
  r_alpha_0.05_sec <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.05, numCores=1)
})

alpha_0.05_par <- system.time({
  r_alpha_0.05_par <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.05, numCores=6)
})

alpha_0.07_sec <- system.time({
  r_alpha_0.07_sec <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.07, numCores=1)
})

alpha_0.07_par <- system.time({
  r_alpha_0.07_par <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.07, numCores=6)
})

alpha_0.1_sec <- system.time({
  r_alpha_0.1_sec <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.1, numCores=1)
})

alpha_0.1_par <- system.time({
  r_alpha_0.1_par <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.1, numCores=6)
})

alpha_0.15_sec <- system.time({
  r_alpha_0.15_sec <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.15, numCores=1)
})

alpha_0.15_par <- system.time({
  r_alpha_0.15_par <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.15, numCores=6)
})

alpha_0.20_sec <- system.time({
  r_alpha_0.20_sec <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.2, numCores=1)
})

alpha_0.20_par <- system.time({
  r_alpha_0.20_par <- getRSKCAlpha(ontMetrics, k=8, L1=3, max_alpha=0.2, numCores=6)
})

cat("Tiempo max_alpha = 0.05 sec:", alpha_0.05_sec["elapsed"], "seg\n")
cat("Tiempo max_alpha = 0.07 sec:", alpha_0.07_sec["elapsed"], "seg\n")
cat("Tiempo max_alpha = 0.1 sec:", alpha_0.1_sec["elapsed"], "seg\n")
cat("Tiempo max_alpha = 0.15 sec:", alpha_0.15_sec["elapsed"], "seg\n")
cat("Tiempo max_alpha = 0.20 sec:", alpha_0.20_sec["elapsed"], "seg\n")

cat("Tiempo max_alpha = 0.05 par:", alpha_0.05_par["elapsed"], "seg\n")
cat("Tiempo max_alpha = 0.07 par:", alpha_0.07_par["elapsed"], "seg\n")
cat("Tiempo max_alpha = 0.1 par:", alpha_0.1_par["elapsed"], "seg\n")
cat("Tiempo max_alpha = 0.15 par:", alpha_0.15_par["elapsed"], "seg\n")
cat("Tiempo max_alpha = 0.20 par:", alpha_0.20_par["elapsed"], "seg\n")

