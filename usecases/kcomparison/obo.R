library(evaluomeR)
source(paste0(wd,"functions.R"))

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
inputPath=paste0(wd,"data/obo-119.csv")
outputDir=paste0(wd,"results-obo")

inputDataObo = read.csv(inputPath, header = TRUE)

plotMetricsMinMax(inputDataObo) # minmax
plotMetricsBoxplot(inputDataObo) # boxplot
plotMetricsCluster(inputDataObo) # cluster
plotMetricsViolin(inputDataObo) # violin

stabilityDataObo <- stabilityRange(data=inputDataObo,
                                   k.range=c(3,15), bs=500,
                                   getImages = FALSE, seed=13606)
qualityDataObo <- qualityRange(data=inputDataObo,
                               k.range=c(3,15), getImages = FALSE,
                               seed=13606)

kOptTableOboFull <- getOptimalKValue(stabilityDataObo, qualityDataObo, k.range=c(3,15))
kOptTableObo <- kOptTableOboFull[, c("Stability_max_k", "Quality_max_k", "Global_optimal_k")]
kOptTableOboValues <- kOptTableOboFull[, c("Stability_max_k_stab", "Quality_max_k_stab")]
csvPath = paste0(outputDir, "/obo.optimal.k",".csv")
#write.csv(kOptTableObo, csvPath, row.names = TRUE)
#View(kOptTableObo)

stabilityTableObo <- getDfStabilityData(stabilityDataObo)
silhouetteTableObo <- getDfQualityData(qualityDataObo)

csvPath = paste0(outputDir, "/obo.stability",".csv")
#write.csv(stabilityTableObo, csvPath, row.names = TRUE)
csvPath = paste0(outputDir, "/obo.goodness",".csv")
#write.csv(silhouetteTableObo, csvPath, row.names = TRUE)

# Conteo de scores estables:

length(which(stabilityTableObo["k_3"]>0.75)) # At least stable, (0.75,0.85]
length(which(silhouetteTableObo["k_3"]>0.5)) # At least reasonable structure, (0.50, 0.70]
