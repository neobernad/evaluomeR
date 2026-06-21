library(evaluomeR)
library(RSKC)
library(sparcl)
library(parallel)

data("ontMetrics")
data("rnaMetrics")
data("bioMetrics")
data("golub")
data("nci60_k8")
data("nci60_k10")

seed=100

# --- Correctness check: serial and parallel must produce identical assay output ---
.stab_correct_sec <- stabilityRange(ontMetrics, k.range=c(3,5), bs=20, numCores=1, seed=100)
.stab_correct_par <- stabilityRange(ontMetrics, k.range=c(3,5), bs=20, numCores=2, seed=100)
stopifnot(identical(
  as.data.frame(assay(.stab_correct_sec)),
  as.data.frame(assay(.stab_correct_par))
))
rm(.stab_correct_sec, .stab_correct_par)
# ---------------------------------------------------------------------------------

#Función auxiliar para series de ejecuciones
ejecutar_experimento_S <- function(dataset, k.range, bs, numCores, seed, nEjec, all_metrics=FALSE){
  tiempos <- numeric(nEjec)
  for(i in 1:nEjec){
    t <- system.time({
      r <- stabilityRange(dataset, k.range=k.range, bs=bs, all_metrics=all_metrics, numCores=numCores, seed=seed)
    })
    tiempos[i] <- t["elapsed"]
  }
  return(list(tiempos = tiempos, media = mean(tiempos)))
}


# CASO 1

# Calculamos el tiempo secuencial
stability1_secuencial <- ejecutar_experimento_S(ontMetrics, k.range=c(3,5), bs=50, numCores=1, seed=100, nEjec=5)

# Calculamos el tiempo paralelo
stability1_parallel_2cores <- ejecutar_experimento_S(ontMetrics, k.range=c(3,5), bs=50, numCores=2, seed=100, nEjec=5)

stability1_parallel_12cores <- ejecutar_experimento_S(ontMetrics, k.range=c(3,5), bs=50, numCores=12, seed=100, nEjec=5)

# Comprobar si los resultados son iguales
# ident1 <- identical(as.data.frame(assay(r_stability1_secuencial[[1]])), as.data.frame(assay(r_stability1_parallel_2cores[[1]])))

cat("Tiempo caso 1 secuencial:", stability1_secuencial$media, "seg\n")
cat("Tiempo caso 1 paralelo 2 cores:", stability1_parallel_2cores$media, "seg\n")
cat("Tiempo caso 1 paralelo 12 cores:", stability1_parallel_12cores$media, "seg\n")


# CASO 2

# Calculamos el tiempo secuencial
stability2_secuencial <- ejecutar_experimento_S(ontMetrics, k.range=c(2,10), bs=100, numCores=1, seed=100, nEjec=5)

# Calculamos el tiempo paralelo
stability2_parallel_2cores <- ejecutar_experimento_S(ontMetrics, k.range=c(2,10), bs=100, numCores=2, seed=100, nEjec=5)

stability2_parallel_12cores <- ejecutar_experimento_S(ontMetrics, k.range=c(2,10), bs=100, numCores=12, seed=100, nEjec=5)

# Comprobar si los resultados son iguales
# ident2 <- identical(as.data.frame(assay(r_stability2_secuencial[[1]])), as.data.frame(assay(r_stability2_parallel_2cores[[1]])))

cat("Tiempo caso 2 secuencial:", stability2_secuencial$media, "seg\n")
cat("Tiempo caso 2 paralelo 2 cores:", stability2_parallel_2cores$media, "seg\n")
cat("Tiempo caso 2 paralelo 12 cores:", stability2_parallel_12cores$media, "seg\n")


# CASO 3

# Calculamos el tiempo secuencial
stability3_secuencial <- ejecutar_experimento_S(golub, k.range=c(2,8), bs=100, numCores=1, seed=100, nEjec=5)

# Calculamos el tiempo paralelo
stability3_parallel_2cores <- ejecutar_experimento_S(golub, k.range=c(2,8), bs=100, numCores=2, seed=100, nEjec=5)

stability3_parallel_12cores <- ejecutar_experimento_S(golub, k.range=c(2,8), bs=100, numCores=12, seed=100, nEjec=5)

# Comprobar si los resultados son iguales
# ident3 <- identical(as.data.frame(assay(r_stability3_secuencial[[1]])), as.data.frame(assay(r_stability3_parallel_2cores[[1]])))

cat("Tiempo caso 3 secuencial:", stability3_secuencial$media, "seg\n")
cat("Tiempo caso 3 paralelo 2 cores:", stability3_parallel_2cores$media, "seg\n")
cat("Tiempo caso 3 paralelo 12 cores:", stability3_parallel_12cores$media, "seg\n")


# CASO 4

# Calculamos el tiempo secuencial
stability4_secuencial <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=1, seed=100, nEjec=5)

# Calculamos el tiempo paralelo
stability4_parallel_2cores <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=2, seed=100, nEjec=5)

stability4_parallel_12cores <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=12, seed=100, nEjec=5)

# Comprobar si los resultados son iguales
# ident4 <- identical(as.data.frame(assay(r_stability4_secuencial[[1]])), as.data.frame(assay(r_stability4_parallel_2cores[[1]])))

cat("Tiempo caso 4 secuencial:", stability4_secuencial$media, "seg\n")
cat("Tiempo caso 4 paralelo 2 cores:", stability4_parallel_2cores$media, "seg\n")
cat("Tiempo caso 4 paralelo 12 cores:", stability4_parallel_12cores$media, "seg\n")


# CASO 5

# Calculamos el tiempo secuencial
stability5_secuencial <- ejecutar_experimento_S(ontMetrics, k.range=c(3,5), bs=50, all_metrics=TRUE, numCores=1, seed=100, nEjec=5)

# Calculamos el tiempo paralelo
stability5_parallel_2cores <- ejecutar_experimento_S(ontMetrics, k.range=c(3,5), bs=50, all_metrics=TRUE, numCores=2, seed=100, nEjec=5)

stability5_parallel_12cores <- ejecutar_experimento_S(ontMetrics, k.range=c(3,5), bs=50, all_metrics=TRUE, numCores=12, seed=100, nEjec=5)

# Comprobar si los resultados son iguales
# ident5 <- identical(as.data.frame(assay(r_stability5_secuencial[[1]])),as.data.frame(assay(r_stability5_parallel_2cores[[1]])))

cat("Tiempo caso 5 secuencial:", stability5_secuencial$media, "seg\n")
cat("Tiempo caso 5 paralelo 2 cores:", stability5_parallel_2cores$media, "seg\n")
cat("Tiempo caso 5 paralelo 12 cores:", stability5_parallel_12cores$media, "seg\n")



# Tests para comprobar parámetro a parámetro la diferencia de rendimiento entre versiones


# NUMERO DE CORES

# Fijamos un dataset (nci60_k8), un número de iteraciones de bootstrapping (bs=100) y rango de k (3,10)

stability_1cores <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=1, seed=100, nEjec=5)

stability_2cores <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=2, seed=100, nEjec=5)

stability_4cores <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=4, seed=100, nEjec=5)

stability_6cores <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=6, seed=100, nEjec=5)

stability_8cores <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=8, seed=100, nEjec=5)

stability_10cores <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=10, seed=100, nEjec=5)

stability_12cores <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=12, seed=100, nEjec=5)

cat("Tiempo secuencial (1 core):", stability_1cores$media, "seg\n")
cat("Tiempo 2 cores:", stability_2cores$media, "seg\n")
cat("Tiempo 4 cores:", stability_4cores$media, "seg\n")
cat("Tiempo 6 cores:", stability_6cores$media, "seg\n")
cat("Tiempo 8 cores:", stability_8cores$media, "seg\n")
cat("Tiempo 10 cores:", stability_10cores$media, "seg\n")
cat("Tiempo 12 cores:", stability_12cores$media, "seg\n")


# NUMERO DE MÉTRICAS

# Fijamos un número de cores (12), un número de iteraciones de bootstrapping (bs=100) y rango de k (3,10)

stability_metrica1_sec <- ejecutar_experimento_S(rnaMetrics, k.range=c(3,10), bs=100, numCores=1, seed=100, nEjec=5)

stability_metrica1_par <- ejecutar_experimento_S(rnaMetrics, k.range=c(3,10), bs=100, numCores=12, seed=100, nEjec=5)

stability_metrica2_sec <- ejecutar_experimento_S(bioMetrics, k.range=c(3,10), bs=100, numCores=1, seed=100, nEjec=5)

stability_metrica2_par <- ejecutar_experimento_S(bioMetrics, k.range=c(3,10), bs=100, numCores=12, seed=100, nEjec=5)

stability_metrica3_sec <- ejecutar_experimento_S(ontMetrics, k.range=c(3,10), bs=100, numCores=1, seed=100, nEjec=5)

stability_metrica3_par <- ejecutar_experimento_S(ontMetrics, k.range=c(3,10), bs=100, numCores=12, seed=100, nEjec=5)

stability_metrica4_sec <- ejecutar_experimento_S(golub, k.range=c(3,10), bs=100, numCores=1, seed=100, nEjec=5)

stability_metrica4_par <- ejecutar_experimento_S(golub, k.range=c(3,10), bs=100, numCores=12, seed=100, nEjec=5)

stability_metrica5_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=1, seed=100, nEjec=5)

stability_metrica5_par <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=12, seed=100, nEjec=5)

stability_metrica6_sec <- ejecutar_experimento_S(nci60_k10, k.range=c(3,10), bs=100, numCores=1, seed=100, nEjec=5)

stability_metrica6_par <- ejecutar_experimento_S(nci60_k10, k.range=c(3,10), bs=100, numCores=12, seed=100, nEjec=5)

cat("Tiempo rnaMetrics sec:", stability_metrica1_sec$media, "seg\n")
cat("Tiempo bioMetrics sec:", stability_metrica2_sec$media, "seg\n")
cat("Tiempo ontMetrics sec:", stability_metrica3_sec$media, "seg\n")
cat("Tiempo golub sec:", stability_metrica4_sec$media, "seg\n")
cat("Tiempo nci60_k8 sec:", stability_metrica5_sec$media, "seg\n")
cat("Tiempo nci60_k10 sec:", stability_metrica6_sec$media, "seg\n")

cat("Tiempo rnaMetrics par:", stability_metrica1_par$media, "seg\n")
cat("Tiempo bioMetrics par:", stability_metrica2_par$media, "seg\n")
cat("Tiempo ontMetrics par:", stability_metrica3_par$media, "seg\n")
cat("Tiempo golub par:", stability_metrica4_par$media, "seg\n")
cat("Tiempo nci60_k8 par:", stability_metrica5_par$media, "seg\n")
cat("Tiempo nci60_k10 par:", stability_metrica6_par$media, "seg\n")


# NÚMERO DE ITERACIONES DE BOOTSTRAPPING

# Fijamos un dataset (nci60_k8), un número de cores (12) y rango de k (3,10)

stability_10bs_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=10, numCores=1, seed=100, nEjec=5)

stability_10bs_par <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=10, numCores=12, seed=100, nEjec=5)

stability_20bs_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=20, numCores=1, seed=100, nEjec=5)

stability_20bs_par <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=20, numCores=12, seed=100, nEjec=5)

stability_50bs_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=50, numCores=1, seed=100, nEjec=5)

stability_50bs_par <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=50, numCores=12, seed=100, nEjec=5)

stability_100bs_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=1, seed=100, nEjec=5)

stability_100bs_par <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=100, numCores=12, seed=100, nEjec=5)

stability_150bs_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=150, numCores=1, seed=100, nEjec=5)

stability_150bs_par <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=150, numCores=12, seed=100, nEjec=5)

stability_200bs_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=200, numCores=1, seed=100, nEjec=5)

stability_200bs_par <- ejecutar_experimento_S(nci60_k8, k.range=c(3,10), bs=200, numCores=12, seed=100, nEjec=5)

cat("Tiempo 10 iteraciones de bootstrapping sec:", stability_10bs_sec$media, "seg\n")
cat("Tiempo 20 iteraciones de bootstrapping sec:", stability_20bs_sec$media, "seg\n")
cat("Tiempo 50 iteraciones de bootstrapping sec:", stability_50bs_sec$media, "seg\n")
cat("Tiempo 100 iteraciones de bootstrapping sec:", stability_100bs_sec$media, "seg\n")
cat("Tiempo 150 iteraciones de bootstrapping sec:", stability_150bs_sec$media, "seg\n")
cat("Tiempo 200 iteraciones de bootstrapping sec:", stability_200bs_sec$media, "seg\n")

cat("Tiempo 10 iteraciones de bootstrapping par:", stability_10bs_par$media, "seg\n")
cat("Tiempo 20 iteraciones de bootstrapping par:", stability_20bs_par$media, "seg\n")
cat("Tiempo 50 iteraciones de bootstrapping par:", stability_50bs_par$media, "seg\n")
cat("Tiempo 100 iteraciones de bootstrapping par:", stability_100bs_par$media, "seg\n")
cat("Tiempo 150 iteraciones de bootstrapping par:", stability_150bs_par$media, "seg\n")
cat("Tiempo 200 iteraciones de bootstrapping par:", stability_200bs_par$media, "seg\n")


# RANGO DE K

# Fijamos un dataset (nci60_k8), un número de cores (12) y número de iteraciones de bootstrapping (bs=100)

stability_k2_3_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(2,3), bs=100, numCores=1, seed=100, nEjec=5)

stability_k2_3_par <- ejecutar_experimento_S(nci60_k8, k.range=c(2,3), bs=100, numCores=12, seed=100, nEjec=5)

stability_k2_5_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(2,5), bs=100, numCores=1, seed=100, nEjec=5)

stability_k2_5_par <- ejecutar_experimento_S(nci60_k8, k.range=c(2,5), bs=100, numCores=12, seed=100, nEjec=5)

stability_k2_8_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(2,8), bs=100, numCores=1, seed=100, nEjec=5)

stability_k2_8_par <- ejecutar_experimento_S(nci60_k8, k.range=c(2,8), bs=100, numCores=12, seed=100, nEjec=5)

stability_k2_10_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(2,10), bs=100, numCores=1, seed=100, nEjec=5)

stability_k2_10_par <- ejecutar_experimento_S(nci60_k8, k.range=c(2,10), bs=100, numCores=12, seed=100, nEjec=5)

stability_k2_12_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(2,12), bs=100, numCores=1, seed=100, nEjec=5)

stability_k2_12_par <- ejecutar_experimento_S(nci60_k8, k.range=c(2,12), bs=100, numCores=12, seed=100, nEjec=5)

stability_k2_15_sec <- ejecutar_experimento_S(nci60_k8, k.range=c(2,15), bs=100, numCores=1, seed=100, nEjec=5)

stability_k2_15_par <- ejecutar_experimento_S(nci60_k8, k.range=c(2,15), bs=100, numCores=12, seed=100, nEjec=5)

cat("Tiempo k = (2,3) sec:", stability_k2_3_sec$media, "seg\n")
cat("Tiempo k = (2,5) sec:", stability_k2_5_sec$media, "seg\n")
cat("Tiempo k = (2,8) sec:", stability_k2_8_sec$media, "seg\n")
cat("Tiempo k = (2,10) sec:", stability_k2_10_sec$media, "seg\n")
cat("Tiempo k = (2,12) sec:", stability_k2_12_sec$media, "seg\n")
cat("Tiempo k = (2,15) sec:", stability_k2_15_sec$media, "seg\n")

cat("Tiempo k = (2,3) par:", stability_k2_3_par$media, "seg\n")
cat("Tiempo k = (2,5) par:", stability_k2_5_par$media, "seg\n")
cat("Tiempo k = (2,8) par:", stability_k2_8_par$media, "seg\n")
cat("Tiempo k = (2,10) par:", stability_k2_10_par$media, "seg\n")
cat("Tiempo k = (2,12) par:", stability_k2_12_par$media, "seg\n")
cat("Tiempo k = (2,15) par:", stability_k2_15_par$media, "seg\n")

