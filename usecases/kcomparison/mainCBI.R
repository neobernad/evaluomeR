library(evaluomeR)
library(DataCombine)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

wd = paste0(dirname(rstudioapi::getSourceEditorContext()$path),"/")
outputDir=paste0(wd,"results-cbi")

evaluomeRSupportedCBI()
optimalKs=NULL
stabValues=NULL
silValues=NULL

cbi_used = c("kmeans","clara","clara_pam")
for (cbi in cbi_used) {
  GLOBAL_CBI = cbi
  optimalKs[[GLOBAL_CBI]] = NULL
  stabValues[[GLOBAL_CBI]] = NULL
  silValues[[GLOBAL_CBI]] = NULL
  source(paste0(wd,"agro.R"))
  source(paste0(wd,"obo.R"))
  optimalKs[[GLOBAL_CBI]]$Agro = as.data.frame(as.numeric(kOptTableAgro$Global_optimal_k),
                                               row.names=rownames(kOptTableObo))
  optimalKs[[GLOBAL_CBI]]$Obo = as.data.frame(as.numeric(kOptTableObo$Global_optimal_k),
                                              row.names=rownames(kOptTableObo))
  colnames(optimalKs[[GLOBAL_CBI]]$Agro) = c("k_opt")
  colnames(optimalKs[[GLOBAL_CBI]]$Obo) = c("k_opt")
  stabValues[[GLOBAL_CBI]]$Agro = stabilityTableAgro
  stabValues[[GLOBAL_CBI]]$Obo = stabilityTableObo
  silValues[[GLOBAL_CBI]]$Agro = silhouetteTableAgro
  silValues[[GLOBAL_CBI]]$Obo = silhouetteTableObo
}

#optimalKs = as.data.frame(optimalKs)
#rownames(optimalKs) = rownames(kOptTableObo)

# Data frame with metrics whose k value matches

#####################################################################################
# Table with stab./qual. scores of optimal k values for each CBI used
#####################################################################################


col_names = append(c("repository", "metric"), cbi_used)

stabilityComparison = setNames(data.frame(matrix(ncol = length(col_names), nrow = 0)), col_names)
qualityComparison = setNames(data.frame(matrix(ncol = length(col_names), nrow = 0)), col_names)
stabCount = 1
silCount = 1

for (metric in rownames(kOptTableObo)) {
  cbi_values = NULL
  for (cbi in cbi_used) { # Get stab. qual. values for each CBI for this metric
    agro_m_kopt = optimalKs[[cbi]]$Agro[metric, "k_opt"]
    obo_m_kopt = optimalKs[[cbi]]$Obo[metric, "k_opt"]
    agro_colname = paste0("k_", agro_m_kopt)
    obo_colname = paste0("k_", obo_m_kopt)
    agro_k_stability = as.numeric(stabValues[[cbi]]$Agro[metric, agro_colname])
    obo_k_stability = as.numeric(stabValues[[cbi]]$Obo[metric, obo_colname])
    agro_k_silhouette = as.numeric(silValues[[cbi]]$Agro[metric, agro_colname])
    obo_k_silhouette = as.numeric(silValues[[cbi]]$Obo[metric, obo_colname])
    cbi_values[[cbi]]$agro_k_stability = agro_k_stability
    cbi_values[[cbi]]$obo_k_stability = obo_k_stability
    cbi_values[[cbi]]$agro_k_silhouette = agro_k_silhouette
    cbi_values[[cbi]]$obo_k_silhouette = obo_k_silhouette
  }

  # Agro Stability
  new_row = c("Agroportal", metric,
              cbi_values[["kmeans"]]$agro_k_stability,
              cbi_values[["clara"]]$agro_k_stability,
              cbi_values[["clara_pam"]]$agro_k_stability
              )
  stabilityComparison <- InsertRow(stabilityComparison, NewRow = new_row, RowNum = stabCount)
  stabCount = stabCount + 1

  # Agro Silhouette
  new_row = c("Agroportal", metric,
              cbi_values[["kmeans"]]$agro_k_silhouette,
              cbi_values[["clara"]]$agro_k_silhouette,
              cbi_values[["clara_pam"]]$agro_k_silhouette              )
  qualityComparison <- InsertRow(qualityComparison, NewRow = new_row, RowNum = silCount)
  silCount = silCount + 1

  # Obo Stability
  new_row = c("OBO Foundry", metric,
              cbi_values[["kmeans"]]$obo_k_stability,
              cbi_values[["clara"]]$obo_k_stability,
              cbi_values[["clara_pam"]]$obo_k_stability
              )
  stabilityComparison <- InsertRow(stabilityComparison, NewRow = new_row, RowNum = stabCount)
  stabCount = stabCount + 1

  # Obo Silhouette
  new_row = c("OBO Foundry", metric,
              cbi_values[["kmeans"]]$obo_k_silhouette,
              cbi_values[["clara"]]$obo_k_silhouette,
              cbi_values[["clara_pam"]]$obo_k_silhouette
              )
  qualityComparison <- InsertRow(qualityComparison, NewRow = new_row, RowNum = silCount)
  silCount = silCount + 1

}

stabilityComparison = na.omit(stabilityComparison)
qualityComparison = na.omit(qualityComparison)

# Column name change from clara_pam to pam
colnames(stabilityComparison)[which(names(stabilityComparison) == "clara_pam")] <- "pam"
colnames(qualityComparison)[which(names(qualityComparison) == "clara_pam")] <- "pam"
measure.vars = c("kmeans", "clara", "pam")


# Stability Agro

agroStabilityComparison = stabilityComparison[stabilityComparison$repository=="Agroportal", ]
agroStabilityComparison$repository <- NULL
agroStabilityComparison = melt(agroStabilityComparison, id="metric", measure.vars = measure.vars)
agroStabilityComparison$value = as.numeric(agroStabilityComparison$value)

agro_sta <- ggplot(data=agroStabilityComparison, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_point(aes(color=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_y_continuous(name="Stability") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.title = element_blank(),
        text = element_text(size=26)
  ) +
  labs(title="Agroportal", subtitle="Stability scores with different algorithms",
       x="Metrics")

# Stability Obo

oboStabilityComparison = stabilityComparison[stabilityComparison$repository=="OBO Foundry", ]
oboStabilityComparison$repository <- NULL
oboStabilityComparison = melt(oboStabilityComparison, id="metric", measure.vars = measure.vars)
oboStabilityComparison$value = as.numeric(oboStabilityComparison$value)

obo_sta <- ggplot(data=oboStabilityComparison, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_point(aes(color=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_y_continuous(name="Stability") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.title = element_blank(),
        text = element_text(size=26)
  ) +
  labs(title="OBO Foundry", subtitle="Stability scores with different algorithms",
       x="Metrics")


# Silhouette Agro

agroSilhouetteComparison = qualityComparison[qualityComparison$repository=="Agroportal", ]
agroSilhouetteComparison$repository <- NULL
agroSilhouetteComparison = melt(agroSilhouetteComparison, id="metric", measure.vars = measure.vars)
agroSilhouetteComparison$value = as.numeric(agroSilhouetteComparison$value)

agro_sil <- ggplot(data=agroSilhouetteComparison, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_point(aes(color=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  scale_y_continuous(name="Goodness") +
  coord_cartesian(ylim = c(0.5, 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.title = element_blank(),
        text = element_text(size=26)
  ) +
  labs(title="Agroportal", subtitle="Goodness scores with different algorithms",
       x="Metrics")

# Silhouette Obo

oboSilhouetteComparison = qualityComparison[qualityComparison$repository=="OBO Foundry", ]
oboSilhouetteComparison$repository <- NULL
oboSilhouetteComparison = melt(oboSilhouetteComparison, id="metric", measure.vars = measure.vars)
oboSilhouetteComparison$value = as.numeric(oboSilhouetteComparison$value)

obo_sil <- ggplot(data=oboSilhouetteComparison, aes(x=metric, y=value, group = variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_point(aes(color=variable)) +
  scale_colour_grey(start = 0, end = 0.5) +
  scale_y_continuous(name="Goodness") + # col="# of cylinders"
  coord_cartesian(ylim = c(0.5, 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        legend.title = element_blank(),
        text = element_text(size=26)
  ) +
  labs(title="OBO Foundry", subtitle="Goodness scores with different algorithms",
       x="Metrics")

# Huge plot!
#grid.arrange(agro_sta, obo_sta, agro_sil, obo_sil, ncol=1)
grid.arrange(agro_sta, agro_sil, ncol=1)
grid.arrange(obo_sta, obo_sil, ncol=1)
