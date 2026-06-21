library(evaluomeR)


data("rnaMetrics")
plotMetricsMinMax(rnaMetrics)
plotMetricsBoxplot(rnaMetrics)
cluster = plotMetricsCluster(ontMetrics, scale = TRUE)
plotMetricsViolin(rnaMetrics)
plotMetricsViolin(ontMetrics, 2)

stabilityData <- stabilityRange(data=rnaMetrics, k.range=c(3,4), bs=20, getImages = FALSE, seed=100)
qualityData <- qualityRange(data=rnaMetrics, k.range=c(3,4), getImages = FALSE, seed=100)

kOptTable <- getOptimalKValue(stabilityData, qualityData, k.range=c(3,4))
kOptTable


df = assay(rnaMetrics)
k.vector1=rep(5,length(colnames(df))-1)
k.vector2=rep(2,length(colnames(df))-1)

plotMetricsClusterComparison(rnaMetrics, k.vector1=k.vector1, k.vector2=k.vector2)
plotMetricsClusterComparison(rnaMetrics, k.vector1=3, k.vector2=c(2,5))
plotMetricsClusterComparison(rnaMetrics, k.vector1=3)

x = as.data.frame(assay(rnaMetrics))

# Multi metric clustering
a = clusterbootWrapper(data=x[c("RIN", "DegFact")], B=100,
                   bootmethod="boot",
                   cbi="kmeans",
                   krange=2, seed=100, gold_standard=NULL)
# bootmean should indicate a stable clustering (Jaccard >= 0.75).
# Exact value is kmeans-RNG-version-specific, so we check the range only.
stopifnot(mean(a$bootmean) >= 0.75 && mean(a$bootmean) <= 1.0)

stab = stability(data=x, k=2, bs=100, seed=100)
stab_mean_RIN <- as.numeric(as.data.frame(assay(stab$stability_mean))[1, "Mean_stability_k_2"])
# stability_mean for RIN should also be in the stable range.
stopifnot(stab_mean_RIN >= 0.75 && stab_mean_RIN <= 1.0)
