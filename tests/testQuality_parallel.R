library(evaluomeR)
library(RSKC)
library(sparcl)
library(parallel)

data("ontMetrics")
data("nci60_k8")

seed=100

# --- Correctness check: serial and parallel must produce identical assay output ---
.qual_correct_sec <- qualityRange(ontMetrics, k.range=c(3,5), quality_index="silhouette", numCores=1, seed=100)
.qual_correct_par <- qualityRange(ontMetrics, k.range=c(3,5), quality_index="silhouette", numCores=2, seed=100)
stopifnot(identical(
  as.data.frame(assay(.qual_correct_sec$k_3)),
  as.data.frame(assay(.qual_correct_par$k_3))
))
stopifnot(identical(
  as.data.frame(assay(.qual_correct_sec$k_5)),
  as.data.frame(assay(.qual_correct_par$k_5))
))
rm(.qual_correct_sec, .qual_correct_par)
# ---------------------------------------------------------------------------------

#Función auxiliar para series de ejecuciones
ejecutar_experimento_Q <- function(dataset, k.range, quality_index, numCores, seed, nEjec, all_metrics=FALSE){
  tiempos <- numeric(nEjec)
  for(i in 1:nEjec){
    t <- system.time({
      r <- qualityRange(dataset, k.range=k.range, all_metrics=all_metrics, quality_index=quality_index, numCores=numCores, seed=seed)
    })
    tiempos[i] <- t["elapsed"]
  }
  return(list(tiempos = tiempos, media = mean(tiempos)))
}


# CASO 1

# Calculamos el tiempo secuencial
quality1_secuencial <- ejecutar_experimento_Q(ontMetrics, k.range=c(3,5), quality_index="silhouette", numCores=1, seed=100, nEjec=5)

# Calculamos el tiempo paralelo
quality1_parallel_2cores <- ejecutar_experimento_Q(ontMetrics, k.range=c(3,5), quality_index="silhouette", numCores=2, seed=100, nEjec=5)

quality1_parallel_12cores <- ejecutar_experimento_Q(ontMetrics, k.range=c(3,5), quality_index="silhouette", numCores=12, seed=100, nEjec=5)

# Comprobar si los resultados son iguales
# ident1 <- identical(as.data.frame(assay(r_quality1_secuencial[[1]])), as.data.frame(assay(r_quality1_parallel_2cores[[1]])))

cat("Tiempo caso 1 secuencial:", quality1_secuencial$media, "seg\n")
cat("Tiempo caso 1 paralelo 2 cores:", quality1_parallel_2cores$media, "seg\n")
cat("Tiempo caso 1 paralelo 12 cores:", quality1_parallel_12cores$media, "seg\n")


# CASO 2

# Calculamos el tiempo secuencial
quality2_secuencial <- ejecutar_experimento_Q(ontMetrics, k.range=c(2,15), quality_index="silhouette", numCores=1, seed=100, nEjec=5)

# Calculamos el tiempo paralelo
quality2_parallel_2cores <- ejecutar_experimento_Q(ontMetrics, k.range=c(2,15), quality_index="silhouette", numCores=2, seed=100, nEjec=5)

quality2_parallel_12cores <- ejecutar_experimento_Q(ontMetrics, k.range=c(2,15), quality_index="silhouette", numCores=12, seed=100, nEjec=5)

# Comprobar si los resultados son iguales
# ident2 <- identical(as.data.frame(assay(r_quality2_secuencial[[1]])), as.data.frame(assay(r_quality2_parallel_2cores[[1]])))

cat("Tiempo caso 2 secuencial:", quality2_secuencial$media, "seg\n")
cat("Tiempo caso 2 paralelo 2 cores:", quality2_parallel_2cores$media, "seg\n")
cat("Tiempo caso 2 paralelo 12 cores:", quality2_parallel_12cores$media, "seg\n")


# CASO 3

# Calculamos el tiempo secuencial
quality3_secuencial <- ejecutar_experimento_Q(nci60_k8, k.range=c(2,8), quality_index="silhouette", numCores=1, seed=100, nEjec=5)

# Calculamos el tiempo paralelo
quality3_parallel_2cores <- ejecutar_experimento_Q(nci60_k8, k.range=c(2,8), quality_index="silhouette", numCores=2, seed=100, nEjec=5)

quality3_parallel_12cores <- ejecutar_experimento_Q(nci60_k8, k.range=c(2,8), quality_index="silhouette", numCores=12, seed=100, nEjec=5)

# Comprobar si los resultados son iguales
# ident3 <- identical(as.data.frame(assay(r_quality3_secuencial[[1]])), as.data.frame(assay(r_quality3_parallel_2cores[[1]])))

cat("Tiempo caso 3 secuencial:", quality3_secuencial$media, "seg\n")
cat("Tiempo caso 3 paralelo 2 cores:", quality3_parallel_2cores$media, "seg\n")
cat("Tiempo caso 3 paralelo 12 cores:", quality3_parallel_12cores$media, "seg\n")


# CASO 4

# Calculamos el tiempo secuencial
quality4_secuencial <- ejecutar_experimento_Q(nci60_k8, k.range=c(3,15), quality_index="silhouette", numCores=1, seed=100, nEjec=5)

# Calculamos el tiempo paralelo
quality4_parallel_2cores <- ejecutar_experimento_Q(nci60_k8, k.range=c(3,15), numCores=2, quality_index="silhouette", seed=100, nEjec=5)

quality4_parallel_12cores <- ejecutar_experimento_Q(nci60_k8, k.range=c(3,15), numCores=12, quality_index="silhouette", seed=100, nEjec=5)

# Comprobar si los resultados son iguales
# ident4 <- identical(as.data.frame(assay(r_quality4_secuencial[[1]])), as.data.frame(assay(r_quality4_parallel_2cores[[1]])))

cat("Tiempo caso 4 secuencial:", quality4_secuencial$media, "seg\n")
cat("Tiempo caso 4 paralelo 2 cores:", quality4_parallel_2cores$media, "seg\n")
cat("Tiempo caso 4 paralelo 12 cores:", quality4_parallel_12cores$media, "seg\n")


# CASO 5

# Calculamos el tiempo secuencial
quality5_secuencial <- ejecutar_experimento_Q(ontMetrics, k.range=c(3,5), all_metrics=TRUE, quality_index="silhouette", numCores=1, seed=100, nEjec=5)

# Calculamos el tiempo paralelo
quality5_parallel_2cores <- ejecutar_experimento_Q(ontMetrics, k.range=c(3,5), all_metrics=TRUE, quality_index="silhouette", numCores=2, seed=100, nEjec=5)

quality5_parallel_12cores <- ejecutar_experimento_Q(ontMetrics, k.range=c(3,5), all_metrics=TRUE, quality_index="silhouette", numCores=12, seed=100, nEjec=5)

# Comprobar si los resultados son iguales
# ident5 <- identical(as.data.frame(assay(r_quality5_secuencial[[1]])), as.data.frame(assay(r_quality5_parallel_2cores[[1]])))

cat("Tiempo caso 5 secuencial:", quality5_secuencial$media, "seg\n")
cat("Tiempo caso 5 paralelo 2 cores:", quality5_parallel_2cores$media, "seg\n")
cat("Tiempo caso 5 paralelo 12 cores:", quality5_parallel_12cores$media, "seg\n")

