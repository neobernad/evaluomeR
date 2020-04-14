library(evaluomeR)

data("rnaMetrics")

dataFrame <- stability(data=rnaMetrics, k=2, bs=100, getImages = FALSE)
dataFrame <- stabilityRange(data=rnaMetrics, k.range=c(2,4), bs=20, getImages = FALSE)
assay(dataFrame)
# > assay(dataFrame)
# Metric    Mean_stability_k_2  Mean_stability_k_3  Mean_stability_k_4
# [1,] "RIN"     "0.825833333333333" "0.778412698412698" "0.69625"
# [2,] "DegFact" "0.955595238095238" "0.977777777777778" "0.820833333333333"
dataFrame <- stabilitySet(data=rnaMetrics, k.set=c(2,3,4), bs=20, getImages = FALSE)

dataFrame <- quality(data=rnaMetrics, k=3, getImages = FALSE)
dataFrame <- qualityRange(data=rnaMetrics, k.range=c(2,15), seed = 20, getImages = FALSE)
assay(getDataQualityRange(dataFrame, 2), 1)
dataFrame1 <- qualitySet(data=rnaMetrics, k.set=c(2,3,4), getImages = FALSE)


dataFrame <- metricsCorrelations(data=rnaMetrics, getImages = FALSE, margins = c(4,4,11,10))
assay(dataFrame, 1)


dataFrame <- stability(data=rnaMetrics, cbi="kmeans", k=2, bs=100, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="clara", k=2, bs=100, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="clara_pam", k=2, bs=100, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="hclust", k=2, bs=100, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="pamk", k=2, bs=100, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="pamk_pam", k=2, bs=100, getImages = FALSE)

dataFrame <- qualityRange(data=rnaMetrics, k.range=c(2,15), getImages = T)
assay(getDataQualityRange(dataFrame, 2), 1)
