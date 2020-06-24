library(evaluomeR)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(stringr)
library(reshape2)
library(ggthemes)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
seed = 13000
plotDir=paste0(wd,"plots")
dataDir=paste0(wd,"results-csv")

dir.create(file.path(plotDir))
dir.create(file.path(dataDir))

getFormattedK <- function(k) {
  return(gsub("^.*_","", k))
}

standardizeStabilityData <- function(stabData, k.range=NULL) {
  stabDf = as.data.frame(assay(stabData)) # Getting first assay, which is 'stabData$stability_mean'
  lengthColnames = length(colnames(stabDf))
  toRemove = list()
  for (i in seq(1, lengthColnames, 1)) {
    colname = colnames(stabDf)[i]
    newColname = gsub("^.*_.*_.*_","k_", colname)
    colnames(stabDf)[i] = newColname
    if (i != 1) { # Skip Metric column
      k = as.numeric(getFormattedK(newColname))
      if (!is.null(k.range) && (k < k.range[1] || k > k.range[2])) {
        toRemove = append(toRemove, newColname)
        next
      }
      stabDf[newColname] = as.numeric(as.character(stabDf[[newColname]]))
    }
  }

  for (columnName in toRemove) {
    stabDf[, columnName] = list(NULL)
    lengthColnames = lengthColnames-1
  }

  inputStartRange = as.numeric(getFormattedK(colnames(stabDf)[2]))
  inputEndRange = as.numeric(getFormattedK(colnames(stabDf)[lengthColnames]))
  if (!is.null(k.range) && (k.range[1] < inputStartRange || k.range[2] > inputEndRange)) {
    # Input k.range is not a subset of the stabData k ranges
    stop("Input k.range [", k.range[1], ", ", k.range[2], "] is not a subset of data range [",
         inputStartRange, ", ", inputEndRange, "]")
  }

  rownames(stabDf) = stabDf$Metric
  stabDf = stabDf[, -1] # Remove "Metric" column, metrics are rownames now
  stabDf <- stabDf[ order(row.names(stabDf)), ]
  return(stabDf)
}

standardizeQualityData <- function(qualData, k.range=NULL) {
  lengthQuality = length(qualData)
  qualRangeStart = getFormattedK(names(qualData)[1])
  qualRangeEnd = getFormattedK(names(qualData)[lengthQuality])
  Metric = NULL
  kValues = list()
  for (i in seq(qualRangeStart, qualRangeEnd, 1)) {
    curQual = as.data.frame(assay(getDataQualityRange(qualData, i)))
    if (i == qualRangeStart) {
      Metric = as.character(curQual$Metric)
    }
    kValues[[i]] = as.numeric(as.character(curQual$Avg_Silhouette_Width))
  }
  qualDf = as.data.frame(Metric)
  for (i in seq(qualRangeStart, qualRangeEnd, 1)) {
    values = kValues[[i]]
    newColname = paste0("k_", i)
    k = as.numeric(getFormattedK(newColname))
    if (!is.null(k.range) && (k < k.range[1] || k > k.range[2])) {
      next
    }
    qualDf[[newColname]] = values
  }

  if (!is.null(k.range) && (k.range[1] < qualRangeStart || k.range[2] > qualRangeEnd)) {
    # Input k.range is not a subset of the stabData k ranges
    stop("Input k.range [", k.range[1], ", ", k.range[2], "] is not a subset of range [",
         qualRangeStart, ", ", qualRangeEnd, "]")
  }

  rownames(qualDf) = qualDf$Metric
  qualDf = qualDf[, -1] # Remove "Metric" column, metrics are rownames now
  qualDf <- qualDf[ order(row.names(qualDf)), ]
  return(qualDf)
}

#### Data loading ----
csai15Path=paste0(wd,"data/data-csai15.csv")
csai16Path=paste0(wd,"data/data-csai16.csv")
csai17Path=paste0(wd,"data/data-csai17.csv")

csis15Path=paste0(wd,"data/data-csis15.csv")
csis16Path=paste0(wd,"data/data-csis16.csv")
csis17Path=paste0(wd,"data/data-csis17.csv")

orms15Path=paste0(wd,"data/data-orms15.csv")
orms16Path=paste0(wd,"data/data-orms16.csv")
orms17Path=paste0(wd,"data/data-orms17.csv")

inputCsai15 = read.csv(csai15Path, header = TRUE, sep = ";")
inputCsai16 = read.csv(csai16Path, header = TRUE, sep = ";")
inputCsai17 = read.csv(csai17Path, header = TRUE, sep = ";")

inputCsis15 = read.csv(csis15Path, header = TRUE, sep = ";")
inputCsis16 = read.csv(csis16Path, header = TRUE, sep = ";")
inputCsis17 = read.csv(csis17Path, header = TRUE, sep = ";")

inputOrms15 = read.csv(orms15Path, header = TRUE, sep = ";")
inputOrms16 = read.csv(orms16Path, header = TRUE, sep = ";")
inputOrms17 = read.csv(orms17Path, header = TRUE, sep = ";")

#### Stability [2,15] ----
k.range=c(2,15)
stabCsai15 <- stabilityRange(data=inputCsai15, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsai16 <- stabilityRange(data=inputCsai16, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsai17 <- stabilityRange(data=inputCsai17, k.range=k.range, bs=100, getImages = FALSE, seed=seed)

stabCsis15 <- stabilityRange(data=inputCsis15, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsis16 <- stabilityRange(data=inputCsis16, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabCsis17 <- stabilityRange(data=inputCsis17, k.range=k.range, bs=100, getImages = FALSE, seed=seed)

stabOrms15 <- stabilityRange(data=inputOrms15, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabOrms16 <- stabilityRange(data=inputOrms16, k.range=k.range, bs=100, getImages = FALSE, seed=seed)
stabOrms17 <- stabilityRange(data=inputOrms17, k.range=k.range, bs=100, getImages = FALSE, seed=seed)

#### Stability SE to DF ----
meanStabCsai15 = standardizeStabilityData(stabCsai15, k.range)
meanStabCsai15$Metric = "CSAI15"
meanStabCsai15$category = "CSAI"
meanStabCsai15$year = "2015"
colnames(meanStabCsai15) = str_remove_all(colnames(meanStabCsai15), "k_")

meanStabCsai16 = standardizeStabilityData(stabCsai16, k.range)
meanStabCsai16$Metric = "CSAI16"
meanStabCsai16$category = "CSAI"
meanStabCsai16$year = "2016"
colnames(meanStabCsai16) = str_remove_all(colnames(meanStabCsai16), "k_")

meanStabCsai17 = standardizeStabilityData(stabCsai17, k.range)
meanStabCsai17$Metric = "CSAI17"
meanStabCsai17$category = "CSAI"
meanStabCsai17$year = "2017"
colnames(meanStabCsai17) = str_remove_all(colnames(meanStabCsai17), "k_")

meanStabCsis15 = standardizeStabilityData(stabCsis15, k.range)
meanStabCsis15$Metric = "CSIS15"
meanStabCsis15$category = "CSIS"
meanStabCsis15$year = "2015"
colnames(meanStabCsis15) = str_remove_all(colnames(meanStabCsis15), "k_")

meanStabCsis16 = standardizeStabilityData(stabCsis16, k.range)
meanStabCsis16$Metric = "CSIS16"
meanStabCsis16$category = "CSIS"
meanStabCsis16$year = "2016"
colnames(meanStabCsis16) = str_remove_all(colnames(meanStabCsis16), "k_")

meanStabCsis17 = standardizeStabilityData(stabCsis17, k.range)
meanStabCsis17$Metric = "CSIS17"
meanStabCsis17$category = "CSIS"
meanStabCsis17$year = "2017"
colnames(meanStabCsis17) = str_remove_all(colnames(meanStabCsis17), "k_")

meanStabOrms15 = standardizeStabilityData(stabOrms15, k.range)
meanStabOrms15$Metric = "ORMS15"
meanStabOrms15$category = "ORMS"
meanStabOrms15$year = "2015"
colnames(meanStabOrms15) = str_remove_all(colnames(meanStabOrms15), "k_")

meanStabOrms16 = standardizeStabilityData(stabOrms16, k.range)
meanStabOrms16$Metric = "ORMS16"
meanStabOrms16$category = "ORMS"
meanStabOrms16$year = "2016"
colnames(meanStabOrms16) = str_remove_all(colnames(meanStabOrms16), "k_")

meanStabOrms17 = standardizeStabilityData(stabOrms17, k.range)
meanStabOrms17$Metric = "ORMS17"
meanStabOrms17$category = "ORMS"
meanStabOrms17$year = "2017"
colnames(meanStabOrms17) = str_remove_all(colnames(meanStabOrms17), "k_")

#### Stability plotting ----

stabDfList = list(meanStabCsai15, meanStabCsai16, meanStabCsai17,
                  meanStabCsis15, meanStabCsis16, meanStabCsis17,
                  meanStabOrms15, meanStabOrms16, meanStabOrms17)

stabPlot = ggplot()
min = Inf
max = -Inf
for (stabDf in stabDfList) {
  stabDf_melt = melt(stabDf, id.vars = c("Metric", "category", "year"))
  colnames(stabDf_melt)[1] = "Journal"
  stabDf_melt$variable = as.integer(stabDf_melt$variable)
  stabPlot = stabPlot +
    geom_line(stabDf_melt,
              mapping = aes(x=variable, y=value, group = 1, linetype=category, colour=category)) +
    geom_point(stabDf_melt,
               mapping = aes(x=variable, y=value, group = 1, shape=year))
  if (min > min(stabDf_melt$value)) {
    min = min(stabDf_melt$value)
  }
  if (max < max(stabDf_melt$value)) {
    max = max(stabDf_melt$value)
  }
}

stabPlot = stabPlot +
  scale_x_continuous(name="k", breaks=1:14, labels=2:15) +
  scale_y_continuous(name="Stability", limits = c(min,max), breaks = c(0.6, 0.75, 0.85), labels=c(0.6, 0.75, 0.85)) +
  scale_colour_grey(start = 0.7, end = 0) +
  theme_calc(base_family = "sans")

ggsave(plot = stabPlot, filename=paste0(plotDir, "/stability_impact_factor.pdf"),
       device="pdf", units="cm", width = 20, height = 10, dpi="print")

#### Quality [2,15] ----
qualCsai15 <- qualityRange(data=inputCsai15, k.range=k.range, getImages = FALSE, seed=seed)
qualCsai16 <- qualityRange(data=inputCsai16, k.range=k.range, getImages = FALSE, seed=seed)
qualCsai17 <- qualityRange(data=inputCsai17, k.range=k.range, getImages = FALSE, seed=seed)

qualCsis15 <- qualityRange(data=inputCsis15, k.range=k.range, getImages = FALSE, seed=seed)
qualCsis16 <- qualityRange(data=inputCsis16, k.range=k.range, getImages = FALSE, seed=seed)
qualCsis17 <- qualityRange(data=inputCsis17, k.range=k.range, getImages = FALSE, seed=seed)

qualOrms15 <- qualityRange(data=inputOrms15, k.range=k.range, getImages = FALSE, seed=seed)
qualOrms16 <- qualityRange(data=inputOrms16, k.range=k.range, getImages = FALSE, seed=seed)
qualOrms17 <- qualityRange(data=inputOrms17, k.range=k.range, getImages = FALSE, seed=seed)

#### Quality SE to DF ----
silCsai15 = standardizeQualityData(qualCsai15, k.range)
silCsai15$Metric = "CSAI15"
silCsai15$category = "CSAI"
silCsai15$year = "2015"
colnames(silCsai15) = str_remove_all(colnames(silCsai15), "k_")

silCsai16 = standardizeQualityData(qualCsai16, k.range)
silCsai16$Metric = "CSAI16"
silCsai16$category = "CSAI"
silCsai16$year = "2016"
colnames(silCsai16) = str_remove_all(colnames(silCsai16), "k_")

silCsai17 = standardizeQualityData(qualCsai17, k.range)
silCsai17$Metric = "CSAI17"
silCsai17$category = "CSAI"
silCsai17$year = "2017"
colnames(silCsai17) = str_remove_all(colnames(silCsai17), "k_")

silCsis15 = standardizeQualityData(qualCsis15, k.range)
silCsis15$Metric = "CSIS15"
silCsis15$category = "CSIS"
silCsis15$year = "2015"
colnames(silCsis15) = str_remove_all(colnames(silCsis15), "k_")

silCsis16 = standardizeQualityData(qualCsis16, k.range)
silCsis16$Metric = "CSIS16"
silCsis16$category = "CSIS"
silCsis16$year = "2016"
colnames(silCsis16) = str_remove_all(colnames(silCsis16), "k_")

silCsis17 = standardizeQualityData(qualCsis17, k.range)
silCsis17$Metric = "CSIS17"
silCsis17$category = "CSIS"
silCsis17$year = "2017"
colnames(silCsis17) = str_remove_all(colnames(silCsis17), "k_")

silOrms15 = standardizeQualityData(qualOrms15, k.range)
silOrms15$Metric = "ORMS15"
silOrms15$category = "ORMS"
silOrms15$year = "2015"
colnames(silOrms15) = str_remove_all(colnames(silOrms15), "k_")

silOrms16 = standardizeQualityData(qualOrms16, k.range)
silOrms16$Metric = "ORMS16"
silOrms16$category = "ORMS"
silOrms16$year = "2016"
colnames(silOrms16) = str_remove_all(colnames(silOrms16), "k_")

silOrms17 = standardizeQualityData(qualOrms17, k.range)
silOrms17$Metric = "ORMS17"
silOrms17$category = "ORMS"
silOrms17$year = "2017"
colnames(silOrms17) = str_remove_all(colnames(silOrms17), "k_")

#### Quality plotting ----

silDfList = list(silCsai15, silCsai16, silCsai17,
                 silCsis15, silCsis16, silCsis17,
                 silOrms15, silOrms16, silOrms17)

silPlot = ggplot()
min = Inf
max = -Inf
for (silDf in silDfList) {
  silDf_melt = melt(silDf, id.vars = c("Metric", "category", "year"))
  colnames(silDf_melt)[1] = "Journal"
  silDf_melt$variable = as.integer(silDf_melt$variable)
  silPlot = silPlot +
    geom_line(silDf_melt,
              mapping = aes(x=variable, y=value, group = 1, linetype=category, colour=category)) +
    geom_point(silDf_melt,
               mapping = aes(x=variable, y=value, group = 1, shape=year))
  if (min > min(silDf_melt$value)) {
    min = min(silDf_melt$value)
  }
  if (max < max(silDf_melt$value)) {
    max = max(silDf_melt$value)
  }
}

silPlot = silPlot +
  scale_x_continuous(name="k", breaks=1:14, labels=2:15) +
  scale_y_continuous(name="Quality", limits = c(min,max), breaks = c(0.25, 0.5, 0.7), labels=c(0.25, 0.5, 0.7)) +
  scale_colour_grey(start = 0.7, end = 0) +
  theme_calc()

ggsave(plot = silPlot, filename=paste0(plotDir, "/silhouette_impact_factor.pdf"),
       device="pdf", units="cm", width = 20, height = 10, dpi="print")

pg = plot_grid(stabPlot, silPlot, align = "v", nrow = 2, rel_heights = c(1/2, 1/2), labels=c("A", "B"))
save_plot(paste0(plotDir, "/stability_silhouette_impact_factor.pdf"), pg, nrow=2, dpi="print")
save_plot(paste0(plotDir, "/stability_silhouette_impact_factor.tiff"), pg, nrow=2, dpi="print", device="tiff")

#### CSV generation ----

# CSAI Stab
write.csv(standardizeStabilityData(stabCsai15, k.range), paste0(dataDir, "/stabilityCsai15.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsai16, k.range), paste0(dataDir, "/stabilityCsai16.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsai17, k.range), paste0(dataDir, "/stabilityCsai17.csv"), row.names = TRUE)

# CSAI Qual
write.csv(standardizeQualityData(qualCsai15, k.range), paste0(dataDir, "/qualityCsai15.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsai16, k.range), paste0(dataDir, "/qualityCsai16.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsai17, k.range), paste0(dataDir, "/qualityCsai17.csv"), row.names = TRUE)

# CSIS Stab
write.csv(standardizeStabilityData(stabCsis15, k.range), paste0(dataDir, "/stabilityCsis15.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsis16, k.range), paste0(dataDir, "/stabilityCsis16.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabCsis17, k.range), paste0(dataDir, "/stabilityCsis17.csv"), row.names = TRUE)

# CSIS Qual
write.csv(standardizeQualityData(qualCsis15, k.range), paste0(dataDir, "/qualityCsis15.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsis16, k.range), paste0(dataDir, "/qualityCsis16.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualCsis17, k.range), paste0(dataDir, "/qualityCsis17.csv"), row.names = TRUE)

# ORMS Stab
write.csv(standardizeStabilityData(stabOrms15, k.range), paste0(dataDir, "/stabilityOrms15.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabOrms16, k.range), paste0(dataDir, "/stabilityOrms16.csv"), row.names = TRUE)
write.csv(standardizeStabilityData(stabOrms17, k.range), paste0(dataDir, "/stabilityOrms17.csv"), row.names = TRUE)

# ORMS Qual
write.csv(standardizeQualityData(qualOrms15, k.range), paste0(dataDir, "/qualityOrms15.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualOrms16, k.range), paste0(dataDir, "/qualityOrms16.csv"), row.names = TRUE)
write.csv(standardizeQualityData(qualOrms17, k.range), paste0(dataDir, "/qualityOrms17.csv"), row.names = TRUE)

#### Optimal k values ----
k.range=c(3,15)

# Csai
kOptCsai15 = getOptimalKValue(stabCsai15, qualCsai15, k.range=k.range)$Global_optimal_k
kOptCsai16 = getOptimalKValue(stabCsai16, qualCsai16, k.range=k.range)$Global_optimal_k
kOptCsai17 = getOptimalKValue(stabCsai17, qualCsai17, k.range=k.range)$Global_optimal_k

# Csis
kOptCsis15 = getOptimalKValue(stabCsis15, qualCsis15, k.range=k.range)$Global_optimal_k
kOptCsis16 = getOptimalKValue(stabCsis16, qualCsis16, k.range=k.range)$Global_optimal_k
kOptCsis17 = getOptimalKValue(stabCsis17, qualCsis17, k.range=k.range)$Global_optimal_k

# Orms
kOptOrms15 = getOptimalKValue(stabOrms15, qualOrms15, k.range=k.range)$Global_optimal_k
kOptOrms16 = getOptimalKValue(stabOrms16, qualOrms16, k.range=k.range)$Global_optimal_k
kOptOrms17 = getOptimalKValue(stabOrms17, qualOrms17, k.range=k.range)$Global_optimal_k

tableKOpt = setNames(as.data.frame(matrix(nrow=3,ncol=3)), c("2015","2016","2017"))
rownames(tableKOpt) <- c("CSAI", "CSIS", "ORMS")
tableKOpt["2015"] <- c(kOptCsai15, kOptCsis15, kOptOrms15)
tableKOpt["2016"] <- c(kOptCsai16, kOptCsis16, kOptOrms16)
tableKOpt["2017"] <- c(kOptCsai17, kOptCsis17, kOptOrms17)
tableKOpt
