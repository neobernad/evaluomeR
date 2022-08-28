library(evaluomeR)
library(RSKC)
library(sparcl)


intputDf = as.data.frame(assay(ontMetrics))
intputDf$Description = NULL

intputDf
annotated_clusters=annotateClustersByMetric(intputDf, k.range=c(2,3), bs=20, seed=100)
plotMetricsCluster(intputDf, k=4)



# due to the ties - there is specific reason to have this be these 3 clusters:
cutree(dend, k = 3)[order.dendrogram(dend)]


inputMatrix = as.matrix(intputDf)
# Note, has to scale matrix otherwise plot (km.perm will not be showing different #non-zero Wjs)
inputMatrix = scale(inputMatrix, TRUE, TRUE)

# Getting tunning paremeter, 1 =< L1 =< sqrt(num.vars)
wbounds = seq(2,sqrt(ncol(inputMatrix)), len=30)
km.perm <- KMeansSparseCluster.permute(x,K=3,wbounds=wbounds,nperms=5)
print(km.perm)
plot(km.perm)
print(km.perm$bestw)


stabilityData <- stabilityRange(data=ontMetrics, k.range=c(3,4), bs=20, getImages = FALSE, seed=100)
qualityData <- qualityRange(data=ontMetrics, k.range=c(3,4), getImages = FALSE, seed=100)
optK <- getOptimalKValue(stabilityData, qualityData, k.range=c(3,4))


for (L1 in c(2,5,10,15)) {
  r3 = RSKC(intputDf, 3, 0.1, L1 = L1, nstart = 200,
            silent=TRUE, scaling = FALSE, correlation = FALSE)
  cat(paste0("L1 value: ", L1,"\n"))
  cat(names(r3$weights)[1], ": ", r3$weights[1],"\n")
  cat(names(r3$weights)[2], ": ", r3$weights[2],"\n")
  cat(names(r3$weights)[3], ": ", r3$weights[3],"\n")
  cat("---\n")
}


data("optd")
truedigit <- rownames(optd)
