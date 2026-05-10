library(evaluomeR)
library(RSKC)
library(sparcl)
library(parallel)

data("ontMetrics")
data("nci60_k8")

seed=100

#Función auxiliar para series de ejecuciones
ejecutar_experimento_Kopt <- function(dataset, k.range, bs, numCores, seed, nEjec, all_metrics=FALSE){
  tiempos <- numeric(nEjec)
  for(i in 1:nEjec){
    t <- system.time({
      rs <- stabilityRange(dataset, k.range=k.range, bs=bs, all_metrics=all_metrics, numCores=numCores, seed=seed)
      rq <- qualityRange(dataset, k.range=k.range, seed=seed)
      k <- getOptimalKValue(rs, rq, k.range=k.range)
    })
    tiempos[i] <- t["elapsed"]
  }
  return(list(tiempos = tiempos, media = mean(tiempos)))
}


# Casos para comprobar la rapidez de getOptimalKValue

# CASO 1

# Calculamos el tiempo de los índices
calculo1_indices <- system.time({
  stabData1 <- stabilityRange(ontMetrics, k.range=c(2,10), bs=100, numCores=12)
  qualData1 <- qualityRange(ontMetrics, k.range=c(2,10))
})

# Calculamos lo que se tarda en realizar el cálculo de k óptimo una vez obtenidos los índices
calculo1_kopt <- system.time({
  kopt1 <- getOptimalKValue(stabData1, qualData1, k.range=c(2,10))
})

cat("Tiempo de cálculo 1 de índices:", calculo1_indices["elapsed"], "seg\n")
cat("Tiempo de cálculo 1 de k óptimo:", calculo1_kopt["elapsed"], "seg\n")


# CASO 2

# Calculamos el tiempo de los índices
calculo2_indices <- system.time({
  stabData2 <- stabilityRange(nci60_k8, k.range=c(2,6), bs=100, numCores=12)
  qualData2 <- qualityRange(nci60_k8, k.range=c(2,6))
})

# Calculamos lo que se tarda en realizar el cálculo de k óptimo una vez obtenidos los índices
calculo2_kopt <- system.time({
  kopt2 <- getOptimalKValue(stabData2, qualData2, k.range=c(2,6))
})

cat("Tiempo de cálculo 2 de índices:", calculo2_indices["elapsed"], "seg\n")
cat("Tiempo de cálculo 2 de k óptimo:", calculo2_kopt["elapsed"], "seg\n")


# CASO 3

# Calculamos el tiempo de los índices
calculo3_indices <- system.time({
  stabData3 <- stabilityRange(nci60_k8, k.range=c(3,10), bs=100, numCores=12)
  qualData3 <- qualityRange(nci60_k8, k.range=c(3,10))
})

# Calculamos lo que se tarda en realizar el cálculo de k óptimo una vez obtenidos los índices
calculo3_kopt <- system.time({
  kopt3 <- getOptimalKValue(stabData3, qualData3, k.range=c(3,10))
})

cat("Tiempo de cálculo 3 de índices:", calculo3_indices["elapsed"], "seg\n")
cat("Tiempo de cálculo 3 de k óptimo:", calculo3_kopt["elapsed"], "seg\n")



# Comparativa cálculo completo k óptimo secuencial y paralelo

# CASO 1

# Calculo secuencial de k óptimo
kopt1_secuencial <- ejecutar_experimento_Kopt(ontMetrics, k.range=c(2,5), bs=50, numCores=1, seed=100, nEjec=5)

# Calculo paralelo de k óptimo
kopt1_paralelo <- ejecutar_experimento_Kopt(ontMetrics, k.range=c(2,5), bs=50, numCores=12, seed=100, nEjec=5)

cat("Tiempo secuencial k_optimo caso 1:", kopt1_secuencial$media, "seg\n")
cat("Tiempo paralelo k_optimo caso 1:", kopt1_paralelo$media, "seg\n")


# CASO 2

# Calculo secuencial de k óptimo
kopt2_secuencial <- ejecutar_experimento_Kopt(ontMetrics, k.range=c(2,10), bs=100, numCores=1, seed=100, nEjec=5)

# Calculo paralelo de k óptimo
kopt2_paralelo <- ejecutar_experimento_Kopt(ontMetrics, k.range=c(2,10), bs=100, numCores=12, seed=100, nEjec=5)

cat("Tiempo secuencial k_optimo caso 2:", kopt2_secuencial$media, "seg\n")
cat("Tiempo paralelo k_optimo caso 2:", kopt2_paralelo$media, "seg\n")


# CASO 3

# Calculo secuencial de k óptimo
kopt3_secuencial <- ejecutar_experimento_Kopt(nci60_k8, k.range=c(2,10), bs=100, numCores=1, seed=100, nEjec=5)

# Calculo paralelo de k óptimo
kopt3_paralelo <- ejecutar_experimento_Kopt(nci60_k8, k.range=c(2,10), bs=100, numCores=12, seed=100, nEjec=5)

cat("Tiempo secuencial k_optimo caso 3:", kopt3_secuencial$media, "seg\n")
cat("Tiempo paralelo k_optimo caso 3:", kopt3_paralelo$media, "seg\n")


# CASO 4

# Calculo secuencial de k óptimo
kopt4_secuencial <- ejecutar_experimento_Kopt(nci60_k8, k.range=c(2,15), bs=150, numCores=1, seed=100, nEjec=5)

# Calculo paralelo de k óptimo
kopt4_paralelo <- ejecutar_experimento_Kopt(nci60_k8, k.range=c(2,15), bs=150, numCores=12, seed=100, nEjec=5)

cat("Tiempo secuencial k_optimo caso 4:", kopt4_secuencial$media, "seg\n")
cat("Tiempo paralelo k_optimo caso 4:", kopt4_paralelo$media, "seg\n")


# Comprobamos si los resultados son iguales
# ident1 <- identical(kopt1_sec, kopt1_par)
# ident2 <- identical(kopt2_sec, kopt2_par)
# ident3 <- identical(kopt3_sec, kopt3_par)
# ident4 <- identical(kopt4_sec, kopt4_par)

# cat("Ident 1:", ident1)
# cat("Ident 2:", ident2)
# cat("Ident 3:", ident3)
# cat("Ident 4:", ident4)

