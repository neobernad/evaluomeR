library(evaluomeR)
library(RSKC)
library(sparcl)
library(parallel)

data("ontMetrics")
data("bioMetrics")
data("rnaMetrics")
data("golub")
data("nci60_k8")

seed=100


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

L1k_2k_ontMetrics <- system.time({
  r_L1k_2k_ontMetrics <- getRSKCL1Boundry(ontMetrics, k=2, clustering="kmeans", seed=100)
})

L1k_4k_ontMetrics <- system.time({
  r_L1k_4k_ontMetrics <- getRSKCL1Boundry(ontMetrics, k=4, clustering="kmeans", seed=100)
})

L1k_8k_ontMetrics <- system.time({
  r_L1k_8k_ontMetrics <- getRSKCL1Boundry(ontMetrics, k=8, clustering="kmeans", seed=100)
})

L1k_15k_ontMetrics <- system.time({
  r_L1k_15k_ontMetrics <- getRSKCL1Boundry(ontMetrics, k=15, clustering="kmeans", seed=100)
})

L1h_ontMetrics <- system.time({
  r_L1h_ontMetrics <- getRSKCL1Boundry(ontMetrics, clustering="hierarchical", seed=100)
})


cat("Tiempo cálculo sparsity bound ontMetrics (k-means clustering, k=2):", L1k_2k_ontMetrics["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound ontMetrics (k-means clustering, k=4):", L1k_4k_ontMetrics["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound ontMetrics (k-means clustering, k=8):", L1k_8k_ontMetrics["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound ontMetrics (k-means clustering, k=15):", L1k_15k_ontMetrics["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound ontMetrics (hierarchical clustering):", L1h_ontMetrics["elapsed"], "seg\n")


cat("Valor sparsity bound ontMetrics (k-means clustering, k=2):", r_L1k_2k_ontMetrics, "\n")
cat("Valor sparsity bound ontMetrics (k-means clustering, k=4):", r_L1k_4k_ontMetrics, "\n")
cat("Valor sparsity bound ontMetrics (k-means clustering, k=8):", r_L1k_8k_ontMetrics, "\n")
cat("Valor sparsity bound ontMetrics (k-means clustering, k=15):", r_L1k_15k_ontMetrics, "\n")
cat("Valor sparsity bound ontMetrics (hierarchical clustering):", r_L1h_ontMetrics, "\n")


# Fijamos dataset = nci60_k8

L1k_2k_nci <- system.time({
  r_L1k_2k_nci <- getRSKCL1Boundry(nci60_k8, k=2, clustering="kmeans", seed=100)
})

L1k_4k_nci <- system.time({
  r_L1k_4k_nci <- getRSKCL1Boundry(nci60_k8, k=4, clustering="kmeans", seed=100)
})

L1k_8k_nci <- system.time({
  r_L1k_8k_nci <- getRSKCL1Boundry(nci60_k8, k=8, clustering="kmeans", seed=100)
})

L1k_15k_nci <- system.time({
  r_L1k_15k_nci <- getRSKCL1Boundry(nci60_k8, k=15, clustering="kmeans", seed=100)
})

L1h_nci <- system.time({
  r_L1h_nci <- getRSKCL1Boundry(nci60_k8, clustering="hierarchical", seed=100)
})


cat("Tiempo cálculo sparsity bound nci60_k8 (k-means clustering, k=2):", L1k_2k_nci["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound nci60_k8 (k-means clustering, k=4):", L1k_4k_nci["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound nci60_k8 (k-means clustering, k=8):", L1k_8k_nci["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound nci60_k8 (k-means clustering, k=15):", L1k_15k_nci["elapsed"], "seg\n")
cat("Tiempo cálculo sparsity bound nci60_k8 (hierarchical clustering):", L1h_nci["elapsed"], "seg\n")


cat("Valor sparsity bound nci60_k8 (k-means clustering, k=2):", r_L1k_2k_nci, "\n")
cat("Valor sparsity bound nci60_k8 (k-means clustering, k=4):", r_L1k_4k_nci, "\n")
cat("Valor sparsity bound nci60_k8 (k-means clustering, k=8):", r_L1k_8k_nci, "\n")
cat("Valor sparsity bound nci60_k8 (k-means clustering, k=15):", r_L1k_15k_nci, "\n")
cat("Valor sparsity bound nci60_k8 (hierarchical clustering):", r_L1h_nci, "\n")
