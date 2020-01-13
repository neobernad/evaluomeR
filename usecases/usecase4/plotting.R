library(evaluomeR)
library(ggplot2)
library(stringr)
library(reshape2)
library(ggthemes)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")

getFormattedK <- function(k) {
  return(gsub("^.*_","", k))
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
matcomp15Path=paste0(wd,"data/data-matcomp15.csv")
matcomp16Path=paste0(wd,"data/data-matcomp16.csv")
matcomp17Path=paste0(wd,"data/data-matcomp17.csv")

multi15Path=paste0(wd,"data/data-multi15.csv")
multi16Path=paste0(wd,"data/data-multi16.csv")
multi17Path=paste0(wd,"data/data-multi17.csv")

multimc15Path=paste0(wd,"data/data-multimc15.csv")
multimc16Path=paste0(wd,"data/data-multimc16.csv")
multimc17Path=paste0(wd,"data/data-multimc17.csv")
outputDir=paste0(wd,"plots")

inputMatcomp15 = read.csv(matcomp15Path, header = TRUE)
inputMatcomp16 = read.csv(matcomp16Path, header = TRUE)
inputMatcomp17 = read.csv(matcomp17Path, header = TRUE)

inputMulti15 = read.csv(multi15Path, header = TRUE)
inputMulti16 = read.csv(multi16Path, header = TRUE)
inputMulti17 = read.csv(multi17Path, header = TRUE)

inputMultimc15 = read.csv(multimc15Path, header = TRUE)
inputMultimc16 = read.csv(multimc16Path, header = TRUE)
inputMultimc17 = read.csv(multimc17Path, header = TRUE)

#### Stability [2,15] ----
stabMatcomp15 <- stabilityRange(data=inputMatcomp15, k.range=c(3,15), bs=100, getImages = FALSE, seed=13606)
stabMatcomp16 <- stabilityRange(data=inputMatcomp16, k.range=c(3,15), bs=100, getImages = FALSE, seed=13606)
stabMatcomp17 <- stabilityRange(data=inputMatcomp17, k.range=c(3,15), bs=100, getImages = FALSE, seed=13606)

stabMulti15 <- stabilityRange(data=inputMulti15, k.range=c(3,15), bs=100, getImages = FALSE, seed=13606)
stabMulti16 <- stabilityRange(data=inputMulti16, k.range=c(3,15), bs=100, getImages = FALSE, seed=13606)
stabMulti17 <- stabilityRange(data=inputMultimc17, k.range=c(3,15), bs=100, getImages = FALSE, seed=13606)

stabMultimc15 <- stabilityRange(data=inputMultimc15, k.range=c(3,15), bs=100, getImages = FALSE, seed=13606)
stabMultimc16 <- stabilityRange(data=inputMultimc16, k.range=c(3,15), bs=100, getImages = FALSE, seed=13606)
stabMultimc17 <- stabilityRange(data=inputMultimc17, k.range=c(3,15), bs=100, getImages = FALSE, seed=13606)

#### Stability SE to DF ----

meanStabMatcomp15 = as.data.frame(assay(stabMatcomp15, "stability_mean"))
meanStabMatcomp15[2:14] <- lapply(meanStabMatcomp15[2:14], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMatcomp15$Metric = "Matcomp15"
meanStabMatcomp15$category = "MCB"
meanStabMatcomp15$year = "2015"
colnames(meanStabMatcomp15) = str_remove_all(colnames(meanStabMatcomp15), "Mean_stability_k_")

meanStabMatcomp16 = as.data.frame(assay(stabMatcomp16, "stability_mean"))
meanStabMatcomp16[2:14] <- lapply(meanStabMatcomp16[2:14], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMatcomp16$Metric = "Matcomp16"
meanStabMatcomp16$category = "MCB"
meanStabMatcomp16$year = "2016"
colnames(meanStabMatcomp16) = str_remove_all(colnames(meanStabMatcomp16), "Mean_stability_k_")

meanStabMatcomp17 = as.data.frame(assay(stabMatcomp17, "stability_mean"))
meanStabMatcomp17[2:14] <- lapply(meanStabMatcomp17[2:14], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMatcomp17$Metric = "Matcomp17"
meanStabMatcomp17$category = "MCB"
meanStabMatcomp17$year = "2017"
colnames(meanStabMatcomp17) = str_remove_all(colnames(meanStabMatcomp17), "Mean_stability_k_")

meanStabMulti15 = as.data.frame(assay(stabMulti15, "stability_mean"))
meanStabMulti15[2:14] <- lapply(meanStabMulti15[2:14], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMulti15$Metric = "Multi15"
meanStabMulti15$category = "MS"
meanStabMulti15$year = "2015"
colnames(meanStabMulti15) = str_remove_all(colnames(meanStabMulti15), "Mean_stability_k_")

meanStabMulti16 = as.data.frame(assay(stabMulti16, "stability_mean"))
meanStabMulti16[2:14] <- lapply(meanStabMulti16[2:14], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMulti16$Metric = "Multi16"
meanStabMulti16$category = "MS"
meanStabMulti16$year = "2016"
colnames(meanStabMulti16) = str_remove_all(colnames(meanStabMulti16), "Mean_stability_k_")

meanStabMulti17 = as.data.frame(assay(stabMulti17, "stability_mean"))
meanStabMulti17[2:14] <- lapply(meanStabMulti17[2:14], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMulti17$Metric = "Multi17"
meanStabMulti17$category = "MS"
meanStabMulti17$year = "2017"
colnames(meanStabMulti17) = str_remove_all(colnames(meanStabMulti17), "Mean_stability_k_")

meanStabMultimc15 = as.data.frame(assay(stabMultimc15, "stability_mean"))
meanStabMultimc15[2:14] <- lapply(meanStabMultimc15[2:14], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMultimc15$Metric = "Multimc15"
meanStabMultimc15$category = "MCB+MS"
meanStabMultimc15$year = "2015"
colnames(meanStabMultimc15) = str_remove_all(colnames(meanStabMultimc15), "Mean_stability_k_")

meanStabMultimc16 = as.data.frame(assay(stabMultimc16, "stability_mean"))
meanStabMultimc16[2:14] <- lapply(meanStabMultimc16[2:14], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMultimc16$Metric = "Multimc16"
meanStabMultimc16$category = "MCB+MS"
meanStabMultimc16$year = "2016"
colnames(meanStabMultimc16) = str_remove_all(colnames(meanStabMultimc16), "Mean_stability_k_")

meanStabMultimc17 = as.data.frame(assay(stabMultimc17, "stability_mean"))
meanStabMultimc17[2:14] <- lapply(meanStabMultimc17[2:14], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMultimc17$Metric = "Multimc17"
meanStabMultimc17$category = "MCB+MS"
meanStabMultimc17$year = "2017"
colnames(meanStabMultimc17) = str_remove_all(colnames(meanStabMultimc17), "Mean_stability_k_")

#### Stability plotting ----

stabDfList = list(meanStabMatcomp15, meanStabMatcomp16, meanStabMatcomp17,
                  meanStabMulti15, meanStabMulti16, meanStabMulti17,
                  meanStabMultimc15, meanStabMultimc16, meanStabMultimc17)

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

stabPlot +
  scale_x_continuous(name="k", breaks=1:13, labels=3:15) +
  scale_y_continuous(name="Stability", limits = c(min,max), breaks = seq(0.5, 1, 0.05), labels=seq(0.5, 1, 0.05)) +
  scale_colour_grey(start = 0.7, end = 0) +
  theme_calc()

#### Quality [3,15] ----
qualMatcomp15 <- qualityRange(data=inputMatcomp15, k.range=c(3,15), getImages = FALSE, seed=13606)
qualMatcomp16 <- qualityRange(data=inputMatcomp16, k.range=c(3,15), getImages = FALSE, seed=13606)
qualMatcomp17 <- qualityRange(data=inputMatcomp17, k.range=c(3,15), getImages = FALSE, seed=13606)

qualMulti15 <- qualityRange(data=inputMulti15, k.range=c(3,15), getImages = FALSE, seed=13606)
qualMulti16 <- qualityRange(data=inputMulti16, k.range=c(3,15), getImages = FALSE, seed=13606)
qualMulti17 <- qualityRange(data=inputMultimc17, k.range=c(3,15), getImages = FALSE, seed=13606)

qualMultimc15 <- qualityRange(data=inputMultimc15, k.range=c(3,15), getImages = FALSE, seed=13606)
qualMultimc16 <- qualityRange(data=inputMultimc16, k.range=c(3,15), getImages = FALSE, seed=13606)
qualMultimc17 <- qualityRange(data=inputMultimc17, k.range=c(3,15), getImages = FALSE, seed=13606)

#### Quality SE to DF ----
k.range=c(3,15)
silMatcomp15 = standardizeQualityData(qualMatcomp15, k.range)
silMatcomp15$Metric = "Matcomp15"
silMatcomp15$category = "MCB"
silMatcomp15$year = "2015"
colnames(silMatcomp15) = str_remove_all(colnames(silMatcomp15), "k_")

silMatcomp16 = standardizeQualityData(qualMatcomp16, k.range)
silMatcomp16$Metric = "Matcomp16"
silMatcomp16$category = "MCB"
silMatcomp16$year = "2016"
colnames(silMatcomp16) = str_remove_all(colnames(silMatcomp16), "k_")

silMatcomp17 = standardizeQualityData(qualMatcomp17, k.range)
silMatcomp17$Metric = "Matcomp17"
silMatcomp17$category = "MCB"
silMatcomp17$year = "2017"
colnames(silMatcomp17) = str_remove_all(colnames(silMatcomp17), "k_")

silMulti15 = standardizeQualityData(qualMulti15, k.range)
silMulti15$Metric = "Multi15"
silMulti15$category = "MS"
silMulti15$year = "2015"
colnames(silMulti15) = str_remove_all(colnames(silMulti15), "k_")

silMulti16 = standardizeQualityData(qualMulti16, k.range)
silMulti16$Metric = "Multi16"
silMulti16$category = "MS"
silMulti16$year = "2016"
colnames(silMulti16) = str_remove_all(colnames(silMulti16), "k_")

silMulti17 = standardizeQualityData(qualMulti17, k.range)
silMulti17$Metric = "Multi17"
silMulti17$category = "MS"
silMulti17$year = "2017"
colnames(silMulti17) = str_remove_all(colnames(silMulti17), "k_")

silMultimc15 = standardizeQualityData(qualMultimc15, k.range)
silMultimc15$Metric = "Multimc15"
silMultimc15$category = "MCB+MS"
silMultimc15$year = "2015"
colnames(silMultimc15) = str_remove_all(colnames(silMultimc15), "k_")

silMultimc16 = standardizeQualityData(qualMultimc16, k.range)
silMultimc16$Metric = "Multimc16"
silMultimc16$category = "MCB+MS"
silMultimc16$year = "2016"
colnames(silMultimc16) = str_remove_all(colnames(silMultimc16), "k_")

silMultimc17 = standardizeQualityData(qualMultimc16, k.range)
silMultimc17$Metric = "Multimc17"
silMultimc17$category = "MCB+MS"
silMultimc17$year = "2017"
colnames(silMultimc17) = str_remove_all(colnames(silMultimc17), "k_")

#### Quality plotting ----

silDfList = list(silMatcomp15, silMatcomp16, silMatcomp17,
                  silMulti15, silMulti16, silMulti17,
                  silMultimc15, silMultimc16, silMultimc17)

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

silPlot +
  scale_x_continuous(name="k", breaks=1:13, labels=3:15) +
  scale_y_continuous(name="Quality", limits = c(min,max), breaks = seq(0.5, 1, 0.05), labels=seq(0.5, 1, 0.05)) +
  scale_colour_grey(start = 0.7, end = 0) +
  theme_calc()

