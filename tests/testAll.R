library(evaluomeR)

data("rnaMetrics")

dataFrame <- stability(data=rnaMetrics, k=4, bs=20, getImages = FALSE)
dataFrame <- stabilityRange(data=rnaMetrics, k.range=c(2,5), bs=20, getImages = FALSE)

dataFrame <- quality(data=rnaMetrics, k=2, getImages = FALSE)
dataFrame <- qualityRange(data=rnaMetrics, k.range=c(2,3), getImages = FALSE)

dataFrame <- correlations(data=rnaMetrics, getImages = FALSE, margins = c(4,4,18,10))
