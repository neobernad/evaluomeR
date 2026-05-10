library(evaluomeR)
library(RSKC)

data("golub")
data("breastCancer")

# golub 

# Gold standard vector

gold_standard_class_G = as.factor(as.vector(unlist(golub["Class"])))

level_mapping <- c("B" = 1, "T" = 2, "M" = 3)
map_strings_to_numbers <- function(strings) {
  factorized <- factor(strings, levels = names(level_mapping))
  return (as.numeric(factorized))
}

gold_standard_vector_G = as.numeric(lapply(gold_standard_class_G, map_strings_to_numbers))

# Limpieza del dataset
r_cleanDataset = evaluomeR::cleanDataset(golub, correlation_threshold = 1)
golub =  r_cleanDataset$dataset

stabData <- stabilityRange(golub, k.range=c(2,6), all_metrics=TRUE, seed=100, gold_standard = gold_standard_vector_G)
# O se puede hacer sin estandarizar con: assay(stabData$cluster_mean)
# Para comprobar los resultados exactos del clustering: assay(stabData$cluster_partition)
rand_G <- standardizeStabilityData(stabData, k.range=c(2,6))
print(rand_G)



# breastCancer

gold_standard_class_B = as.factor(as.vector(unlist(breastCancer["diagnosis"])))

level_mapping <- c("B" = 1, "M" = 2)
map_strings_to_numbers <- function(strings) {
  factorized <- factor(strings, levels = names(level_mapping))
  return (as.numeric(factorized))
}

gold_standard_vector_B = as.numeric(lapply(gold_standard_class_B, map_strings_to_numbers))

# Limpieza del dataset
breastCancer_clean <- breastCancer[, -2]

stabData <- stabilityRange(breastCancer_clean, k.range=c(2,6), all_metrics=TRUE, seed=100, gold_standard = gold_standard_vector_B)
# O se puede hacer sin estandarizar con: assay(stabData$cluster_mean)
# Para comprobar los resultados exactos del clustering: assay(stabData$cluster_partition)
rand_B <- standardizeStabilityData(stabData, k.range=c(2,6))
print(rand_B)

