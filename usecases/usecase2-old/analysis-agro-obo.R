library(evaluomeR)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
inputAgroPath=paste0(wd,"data/agro.csv")
inputOboPath=paste0(wd,"data/obo-119.csv")
outputDir=paste0(wd,"results-agro-obo")

agroData = read.csv(inputAgroPath, header = TRUE)
oboData = read.csv(inputOboPath, header = TRUE)
inputCSV = rbind(agroData, oboData)

# Stability [2,6]
dataFrame <- stabilityRange(data=inputData, k.range=c(2,6), bs=100, getImages = TRUE)
dataStability = assay(dataFrame, 1)
csvPath = paste0(outputDir, "/stability.csv")
write.csv(dataStability, csvPath, row.names = FALSE)

# Quality [2,6]
dataFrame <- qualityRange(data=inputData, k.range=c(2,6), getImages = TRUE)
for (i in 2:6) {
  data = assay(getDataQualityRange(dataFrame, i), 1)
  csvPath = paste0(outputDir, "/quality_k=", i,".csv")
  write.csv(data, csvPath, row.names = FALSE)
}

# Correlations
dataFrame <- metricsCorrelations(inputData, margins=c(0,6,6,5), getImages=TRUE)
csvPath = paste0(outputDir, "/agro.obo.correlation",".csv")
write.csv(assay(dataFrame, 1), csvPath, row.names = FALSE)
