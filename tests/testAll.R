library(evaluomeR)

rnaMetrics <- loadSample("rna-metrics")
location <- NULL # "/home/user/outputdir"


dataFrame <- stability(data=rnaMetrics, k=4, bs=20, getImages = FALSE, path=location)
dataFrame <- stabilityRange(data=rnaMetrics, k.range=c(2,5), bs=20, getImages = FALSE, path=location)

dataFrame <- quality(data=rnaMetrics, k=2, label="Exp.1:", getImages = FALSE, path=location)

dataFrame <- qualityRange(data=rnaMetrics, label="File1:", k.range=c(2,3), getImages = FALSE, path=location)
dataFrame <- correlations(data=rnaMetrics, getImages = FALSE, margins = c(4,4,12,10), path=location)
