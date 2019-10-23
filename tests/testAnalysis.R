library(evaluomeR)

data("rnaMetrics")
plotMetricsMinMax(rnaMetrics)
plotMetricsBoxplot(rnaMetrics)
plotMetricsCluster(rnaMetrics)

set.seed(100)
stabilityData <- stabilityRange(data=rnaMetrics, k.range=c(2,4), bs=20, getImages = FALSE)
set.seed(100)
qualityData <- qualityRange(data=rnaMetrics, k.range=c(2,4), getImages = FALSE)

kOptTable <- getOptimalKValue(stabilityData, qualityData, k.range=c(2,4))
kOptTable
