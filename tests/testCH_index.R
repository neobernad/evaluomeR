library(evaluomeR)
library(RSKC)
library(sparcl)
library(fpc)

data("ontMetrics")
data("nci60_k8")

seed=100

#Función auxiliar para series de ejecuciones
ejecutar_experimento_Quality_indices <- function(dataset, k, numCores, seed, quality_index, nEjec, all_metrics=FALSE){
  tiempos <- numeric(nEjec)
  resultados_calidad <- numeric(nEjec)
  for(i in 1:nEjec){
    t <- system.time({
      r <- quality(dataset, k=k, all_metrics=all_metrics, quality_index=quality_index, numCores=numCores, seed=seed)
    })
    if(quality_index=="ch"){
      resultados_calidad[i] = mean(r$CH)
    }else{
      if(numCores == 1){
        qual = evaluomeR::standardizeQualityData(r, k.range = c(k,k))
        resultados_calidad[i] = mean(qual[[1]])
      }else{
        resultados_calidad[i] = mean(r$Silhouette)
      }
    }
    tiempos[i] <- t["elapsed"]
  }
  return(list(tiempos = tiempos, media_tiempos = mean(tiempos), media_calidad = mean(resultados_calidad)))
}


# ontMetrics

quality1_2k_silhouette <- ejecutar_experimento_Quality_indices(ontMetrics, k=2, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality1_2k_calinhara <- ejecutar_experimento_Quality_indices(ontMetrics, k=2, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality1_4k_silhouette <- ejecutar_experimento_Quality_indices(ontMetrics, k=4, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality1_4k_calinhara <- ejecutar_experimento_Quality_indices(ontMetrics, k=4, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality1_6k_silhouette <- ejecutar_experimento_Quality_indices(ontMetrics, k=6, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality1_6k_calinhara <- ejecutar_experimento_Quality_indices(ontMetrics, k=6, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality1_8k_silhouette <- ejecutar_experimento_Quality_indices(ontMetrics, k=8, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality1_8k_calinhara <- ejecutar_experimento_Quality_indices(ontMetrics, k=8, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality1_10k_silhouette <- ejecutar_experimento_Quality_indices(ontMetrics, k=10, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality1_10k_calinhara <- ejecutar_experimento_Quality_indices(ontMetrics, k=10, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality1_12k_silhouette <- ejecutar_experimento_Quality_indices(ontMetrics, k=12, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality1_12k_calinhara <- ejecutar_experimento_Quality_indices(ontMetrics, k=12, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality1_15k_silhouette <- ejecutar_experimento_Quality_indices(ontMetrics, k=15, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality1_15k_calinhara <- ejecutar_experimento_Quality_indices(ontMetrics, k=15, numCores=1, seed=100, quality_index = "ch", nEjec=5)

cat("Tiempo ejecución cálculo de Silhouette en ontMetrics k = 2:", quality1_2k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Calinski-Harabasz en ontMetrics k = 2:", quality1_2k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Silhouette en ontMetrics k = 4:", quality1_4k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Calinski-Harabasz en ontMetrics k = 4:", quality1_4k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Silhouette en ontMetrics k = 6:", quality1_6k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Calinski-Harabasz en ontMetrics k = 6:", quality1_6k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Silhouette en ontMetrics k = 8:", quality1_8k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Calinski-Harabasz en ontMetrics k = 8:", quality1_8k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Silhouette en ontMetrics k = 10:", quality1_10k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Calinski-Harabasz en ontMetrics k = 10:", quality1_10k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Silhouette en ontMetrics k = 12:", quality1_12k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Calinski-Harabasz en ontMetrics k = 12:", quality1_12k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Silhouette en ontMetrics k = 15:", quality1_15k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución cálculo de Calinski-Harabasz en ontMetrics k = 15:", quality1_15k_calinhara$media_tiempos, "\n")

cat("Calidad cálculo de Silhouette en ontMetrics k = 2:", quality1_2k_silhouette$media_calidad, "\n")
cat("Calidad cálculo de Calinski-Harabasz en ontMetrics k = 2:", quality1_2k_calinhara$media_calidad, "\n")
cat("Calidad cálculo de Silhouette en ontMetrics k = 4:", quality1_4k_silhouette$media_calidad, "\n")
cat("Calidad cálculo de Calinski-Harabasz en ontMetrics k = 4:", quality1_4k_calinhara$media_calidad, "\n")
cat("Calidad cálculo de Silhouette en ontMetrics k = 6:", quality1_6k_silhouette$media_calidad, "\n")
cat("Calidad cálculo de Calinski-Harabasz en ontMetrics k = 6:", quality1_6k_calinhara$media_calidad, "\n")
cat("Calidad cálculo de Silhouette en ontMetrics k = 8:", quality1_8k_silhouette$media_calidad, "\n")
cat("Calidad cálculo de Calinski-Harabasz en ontMetrics k = 8:", quality1_8k_calinhara$media_calidad, "\n")
cat("Calidad cálculo de Silhouette en ontMetrics k = 10:", quality1_10k_silhouette$media_calidad, "\n")
cat("Calidad cálculo de Calinski-Harabasz en ontMetrics k = 10:", quality1_10k_calinhara$media_calidad, "\n")
cat("Calidad cálculo de Silhouette en ontMetrics k = 12:", quality1_12k_silhouette$media_calidad, "\n")
cat("Calidad cálculo de Calinski-Harabasz en ontMetrics k = 12:", quality1_12k_calinhara$media_calidad, "\n")
cat("Calidad cálculo de Silhouette en ontMetrics k = 15:", quality1_15k_silhouette$media_calidad, "\n")
cat("Calidad cálculo de Calinski-Harabasz en ontMetrics k = 15:", quality1_15k_calinhara$media_calidad, "\n")


# nci60_k8

quality2_2k_silhouette <- ejecutar_experimento_Quality_indices(nci60_k8, k=2, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality2_2k_calinhara <- ejecutar_experimento_Quality_indices(nci60_k8, k=2, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality2_4k_silhouette <- ejecutar_experimento_Quality_indices(nci60_k8, k=4, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality2_4k_calinhara <- ejecutar_experimento_Quality_indices(nci60_k8, k=4, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality2_6k_silhouette <- ejecutar_experimento_Quality_indices(nci60_k8, k=6, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality2_6k_calinhara <- ejecutar_experimento_Quality_indices(nci60_k8, k=6, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality2_8k_silhouette <- ejecutar_experimento_Quality_indices(nci60_k8, k=8, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality2_8k_calinhara <- ejecutar_experimento_Quality_indices(nci60_k8, k=8, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality2_10k_silhouette <- ejecutar_experimento_Quality_indices(nci60_k8, k=10, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality2_10k_calinhara <- ejecutar_experimento_Quality_indices(nci60_k8, k=10, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality2_12k_silhouette <- ejecutar_experimento_Quality_indices(nci60_k8, k=12, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality2_12k_calinhara <- ejecutar_experimento_Quality_indices(nci60_k8, k=12, numCores=1, seed=100, quality_index = "ch", nEjec=5)

quality2_15k_silhouette <- ejecutar_experimento_Quality_indices(nci60_k8, k=15, numCores=1, seed=100, quality_index = "silhouette", nEjec=5)

quality2_15k_calinhara <- ejecutar_experimento_Quality_indices(nci60_k8, k=15, numCores=1, seed=100, quality_index = "ch", nEjec=5)

cat("Tiempo ejecución índice de Silhouette en nci60_k8 k = 2:", quality2_2k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución índice de Calinski-Harabasz en nci60_k8 k = 2:", quality2_2k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución índice de Silhouette en nci60_k8 k = 4:", quality2_4k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución índice de Calinski-Harabasz en nci60_k8 k = 4:", quality2_4k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución índice de Silhouette en nci60_k8 k = 6:", quality2_6k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución índice de Calinski-Harabasz en nci60_k8 k = 6:", quality2_6k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución índice de Silhouette en nci60_k8 k = 8:", quality2_8k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución índice de Calinski-Harabasz en nci60_k8 k = 8:", quality2_8k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución índice de Silhouette en nci60_k8 k = 10:", quality2_10k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución índice de Calinski-Harabasz en nci60_k8 k = 10:", quality2_10k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución índice de Silhouette en nci60_k8 k = 12:", quality2_12k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución índice de Calinski-Harabasz en nci60_k8 k = 12:", quality2_12k_calinhara$media_tiempos, "\n")
cat("Tiempo ejecución índice de Silhouette en nci60_k8 k = 15:", quality2_15k_silhouette$media_tiempos, "\n")
cat("Tiempo ejecución índice de Calinski-Harabasz en nci60_k8 k = 15:", quality2_15k_calinhara$media_tiempos, "\n")

cat("Calidad índice de Silhouette en nci60_k8 k = 2:", quality2_2k_silhouette$media_calidad, "\n")
cat("Calidad índice de Calinski-Harabasz en nci60_k8 k = 2:", quality2_2k_calinhara$media_calidad, "\n")
cat("Calidad índice de Silhouette en nci60_k8 k = 4:", quality2_4k_silhouette$media_calidad, "\n")
cat("Calidad índice de Calinski-Harabasz en nci60_k8 k = 4:", quality2_4k_calinhara$media_calidad, "\n")
cat("Calidad índice de Silhouette en nci60_k8 k = 6:", quality2_6k_silhouette$media_calidad, "\n")
cat("Calidad índice de Calinski-Harabasz en nci60_k8 k = 6:", quality2_6k_calinhara$media_calidad, "\n")
cat("Calidad índice de Silhouette en nci60_k8 k = 8:", quality2_8k_silhouette$media_calidad, "\n")
cat("Calidad índice de Calinski-Harabasz en nci60_k8 k = 8:", quality2_8k_calinhara$media_calidad, "\n")
cat("Calidad índice de Silhouette en nci60_k8 k = 10:", quality2_10k_silhouette$media_calidad, "\n")
cat("Calidad índice de Calinski-Harabasz en nci60_k8 k = 10:", quality2_10k_calinhara$media_calidad, "\n")
cat("Calidad índice de Silhouette en nci60_k8 k = 12:", quality2_12k_silhouette$media_calidad, "\n")
cat("Calidad índice de Calinski-Harabasz en nci60_k8 k = 12:", quality2_12k_calinhara$media_calidad, "\n")
cat("Calidad índice de Silhouette en nci60_k8 k = 15:", quality2_15k_silhouette$media_calidad, "\n")
cat("Calidad índice de Calinski-Harabasz en nci60_k8 k = 15:", quality2_15k_calinhara$media_calidad, "\n")

