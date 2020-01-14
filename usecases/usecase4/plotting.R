library(evaluomeR)
library(ggplot2)
library(stringr)
library(reshape2)
library(ggthemes)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")

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
inputAgroPath=paste0(wd,"data/agro.csv")
inputOboPath=paste0(wd,"data/obo-119.csv")
outputDir=paste0(wd,"plots")

agroData = read.csv(inputAgroPath, header = TRUE)
oboData = read.csv(inputOboPath, header = TRUE)
bothData = rbind(agroData, oboData)


#### Stability [2,6] ----
stabAgro <- stabilityRange(data=agroData, k.range=c(2,6), bs=100, getImages = FALSE, seed=13606)
stabObo <- stabilityRange(data=oboData, k.range=c(2,6), bs=100, getImages = FALSE, seed=13606)
stabBoth <- stabilityRange(data=bothData, k.range=c(2,6), bs=100, getImages = FALSE, seed=13606)

#### Stability SE to DF ----
k.range=c(2,6)
meanStabAgro = standardizeStabilityData(stabAgro, k.range)
meanStabAgro$Metric = rownames(meanStabAgro)
meanStabAgro$Repository = "AGRO"
colnames(meanStabAgro) = str_remove_all(colnames(meanStabAgro), "k_")

meanStabObo = standardizeStabilityData(stabObo, k.range)
meanStabObo$Metric = rownames(meanStabObo)
meanStabObo$Repository = "OBO"
colnames(meanStabObo) = str_remove_all(colnames(meanStabObo), "k_")

meanStabBoth = standardizeStabilityData(stabBoth, k.range)
meanStabBoth$Metric = rownames(meanStabBoth)
meanStabBoth$Repository = "AGRO+OBO"
colnames(meanStabBoth) = str_remove_all(colnames(meanStabBoth), "k_")

#### Stability plotting ----

stabDfList = list(meanStabAgro, meanStabObo, meanStabBoth)

stabPlot = ggplot()
min = Inf
max = -Inf
for (stabDf in stabDfList) {
  stabDf = stabDf[which(stabDf$Metric == "ANOnto" | stabDf$Metric == "CBOOnto"), ]
  stabDf_melt = melt(stabDf, id.vars = c("Metric", "Repository"))
  stabDf_melt$variable = as.integer(stabDf_melt$variable)
  for (metric in c("ANOnto", "CBOOnto")) {
    current_melt = stabDf_melt[which(stabDf_melt$Metric == metric), ]
    stabPlot = stabPlot +
      geom_line(current_melt,
                mapping = aes(x=variable, y=value, group = 1, linetype=Metric, colour=Metric)) +
      geom_point(stabDf_melt,
                mapping = aes(x=variable, y=value, group = 1, shape=Repository))
  }

  if (min > min(stabDf_melt$value)) {
    min = min(stabDf_melt$value)
  }
  if (max < max(stabDf_melt$value)) {
    max = max(stabDf_melt$value)
  }
}

stabPlot +
  scale_x_continuous(name="k", breaks=1:5, labels=2:6) +
  scale_y_continuous(name="Stability", limits = c(min,max), breaks = seq(0.5, 1, 0.025), labels=seq(0.5, 1, 0.025)) +
  scale_colour_grey(start = 0.7, end = 0) +
  theme_calc()


