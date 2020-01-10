library(evaluomeR)
library(ggplot2)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
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
stabMatcomp15 <- stabilityRange(data=inputMatcomp15, k.range=c(2,15), bs=100, getImages = TRUE, seed=13606)
stabMatcomp16 <- stabilityRange(data=inputMatcomp16, k.range=c(2,15), bs=100, getImages = TRUE, seed=13606)
stabMatcomp17 <- stabilityRange(data=inputMatcomp17, k.range=c(2,15), bs=100, getImages = TRUE, seed=13606)

stabMulti15 <- stabilityRange(data=inputMulti15, k.range=c(2,15), bs=100, getImages = TRUE, seed=13606)
stabMulti16 <- stabilityRange(data=inputMulti16, k.range=c(2,15), bs=100, getImages = TRUE, seed=13606)
stabMulti17 <- stabilityRange(data=inputMultimc17, k.range=c(2,15), bs=100, getImages = TRUE, seed=13606)

stabMultimc15 <- stabilityRange(data=inputMultimc15, k.range=c(2,15), bs=100, getImages = TRUE, seed=13606)
stabMultimc16 <- stabilityRange(data=inputMultimc16, k.range=c(2,15), bs=100, getImages = TRUE, seed=13606)
stabMultimc17 <- stabilityRange(data=inputMultimc17, k.range=c(2,15), bs=100, getImages = TRUE, seed=13606)

#### Stability SE to DF ----

meanStabMatcomp15 = as.data.frame(assay(stabMatcomp15, "stability_mean"))
meanStabMatcomp15[2:15] <- lapply(meanStabMatcomp15[2:15], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMatcomp15$Metric = "Matcomp15"
colnames(meanStabMatcomp15) = str_remove_all(colnames(meanStabMatcomp15), "Mean_stability_k_")

meanStabMatcomp16 = as.data.frame(assay(stabMatcomp16, "stability_mean"))
meanStabMatcomp16[2:15] <- lapply(meanStabMatcomp16[2:15], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMatcomp16$Metric = "Matcomp16"
colnames(meanStabMatcomp16) = str_remove_all(colnames(meanStabMatcomp16), "Mean_stability_k_")

meanStabMatcomp17 = as.data.frame(assay(stabMatcomp17, "stability_mean"))
meanStabMatcomp17[2:15] <- lapply(meanStabMatcomp17[2:15], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMatcomp17$Metric = "Matcomp17"
colnames(meanStabMatcomp17) = str_remove_all(colnames(meanStabMatcomp17), "Mean_stability_k_")

meanStabMulti15 = as.data.frame(assay(stabMulti15, "stability_mean"))
meanStabMulti15[2:15] <- lapply(meanStabMulti15[2:15], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMulti15$Metric = "Multi15"
colnames(meanStabMulti15) = str_remove_all(colnames(meanStabMulti15), "Mean_stability_k_")

meanStabMulti16 = as.data.frame(assay(stabMulti16, "stability_mean"))
meanStabMulti16[2:15] <- lapply(meanStabMulti16[2:15], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMulti16$Metric = "Multi16"
colnames(meanStabMulti16) = str_remove_all(colnames(meanStabMulti16), "Mean_stability_k_")

meanStabMulti17 = as.data.frame(assay(stabMulti17, "stability_mean"))
meanStabMulti17[2:15] <- lapply(meanStabMulti17[2:15], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMulti17$Metric = "Multi17"
colnames(meanStabMulti17) = str_remove_all(colnames(meanStabMulti17), "Mean_stability_k_")

meanStabMultimc15 = as.data.frame(assay(stabMultimc15, "stability_mean"))
meanStabMultimc15[2:15] <- lapply(meanStabMultimc15[2:15], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMultimc15$Metric = "Multimc15"
colnames(meanStabMultimc15) = str_remove_all(colnames(meanStabMultimc15), "Mean_stability_k_")

meanStabMultimc16 = as.data.frame(assay(stabMultimc16, "stability_mean"))
meanStabMultimc16[2:15] <- lapply(meanStabMultimc16[2:15], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMultimc16$Metric = "Multimc16"
colnames(meanStabMultimc16) = str_remove_all(colnames(meanStabMultimc16), "Mean_stability_k_")

meanStabMultimc17 = as.data.frame(assay(stabMultimc17, "stability_mean"))
meanStabMultimc17[2:15] <- lapply(meanStabMultimc17[2:15], function(x) {if(is.factor(x)) as.numeric(as.character(x)) else x})
meanStabMultimc17$Metric = "Multimc17"
colnames(meanStabMultimc17) = str_remove_all(colnames(meanStabMultimc17), "Mean_stability_k_")

#### Stability plotting ----

stabDfList = list(meanStabMatcomp15, meanStabMatcomp16)

stabPlot = ggplot()

for (stabDf in stabDfList) {
  stabDf_melt = melt(stabDf, id.vars = "Metric")
  colnames(stabDf_melt)[1] = "Journal"
  stabDf_melt$variable = as.integer(stabDf_melt$variable)
  stabPlot = stabPlot +
    geom_line(stabDf_melt,
              mapping = aes(x=variable, y=value, group = 1, colour=Journal)) +
    geom_point(stabDf_melt,
               mapping = aes(x=variable, y=value, group = 1, colour=Journal))
}

stabPlot +
  scale_x_continuous(name="k", breaks=1:14, labels=2:15) +
  ylim(0,1)
