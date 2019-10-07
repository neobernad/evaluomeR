library(evaluomeR)

data("rnaMetrics")

dataFrame <- stability(data=rnaMetrics, k=2, bs=20, getImages = FALSE)
dataFrame <- stabilityRange(data=rnaMetrics, k.range=c(2,3), bs=20, getImages = FALSE)

dataFrame <- quality(data=rnaMetrics, k=3, getImages = T)
dataFrame <- qualityRange(data=rnaMetrics, k.range=c(2,3), getImages = T)
assay(getDataQualityRange(dataFrame, 2), 1)

dataFrame <- metricsCorrelations(data=rnaMetrics, getImages = FALSE, margins = c(4,4,11,10))
assay(dataFrame, 1)
