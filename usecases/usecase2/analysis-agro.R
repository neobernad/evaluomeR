library(evaluomeR)

wd = paste0(getwd(),"/")
inputPath=paste0(wd,"data/agro.csv")
outputDir=paste0(wd,"results-agro")

inputData = read.csv(inputPath, header = TRUE)

dataFrame <- stabilityRange(data=inputData, k.range=c(2,6), bs=20, getImages = T)
dataFrame <- qualityRange(data=inputData, k.range=c(2,6), getImages = T)
dataStability = assay(dataFrame, 1)
