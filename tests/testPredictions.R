library(evaluomeR)

globalMetricDf = globalMetric(ontMetrics, k.range = c(2,3), nrep=10, criterion="AIC", PCA=F)
globalMetricDf

