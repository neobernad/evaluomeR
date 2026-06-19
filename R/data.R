#' @title Dataset: Structural ontology metrics
#' @description
#' Structural ontology metrics, 19 metrics measuring
#' structural aspects of bio-ontologies have been analysed on two different
#' corpora of ontologies: OBO Foundry and AgroPortal
#'
#' @references
#' \insertRef{ontoeval}{evaluomeR}
#' @usage data("ontMetrics")
"ontMetrics"

#' @title Dataset: RNA quality metrics
#' @description
#' RNA quality metrics for the assessment of gene expression
#' differences, 2 quality metrics from 16 aliquots of a unique batch of RNA
#' Samples. The metrics are Degradation Factor (DegFact) and RNA Integrity Number
#' (RIN)
#'
#' @references
#' \insertRef{imbeaud2005towards}{evaluomeR}
#' @usage data("rnaMetrics")
"rnaMetrics"

#' @title Dataset: Metrics for biological pathways
#' @description
#' Metrics for biological pathways, 2 metrics that
#' quantitative characterizations of the importance of regulation in biochemical
#' pathway systems, including systems designed for applications in synthetic
#' biology or metabolic engineering. The metrics are reachability and efficiency
#'
#' @references
#' \insertRef{davis2018metrics}{evaluomeR}
#' @usage data("bioMetrics")
"bioMetrics"

#' @title Dataset: Structural ontology metrics (OBO Foundry subset)
#' @description
#' Structural ontology metrics measured on the OBO Foundry corpus of
#' bio-ontologies. A subset of the ontMetrics dataset restricted to OBO
#' Foundry ontologies.
#' @usage data("ontMetricsOBO")
"ontMetricsOBO"

#' @title Dataset: NCI-60 cancer cell line metrics (k=8 partition)
#' @description
#' Bioinformatics metrics computed on the NCI-60 cancer cell line panel,
#' pre-partitioned into 8 clusters. Used for benchmarking clustering methods.
#' @usage data("nci60_k8")
"nci60_k8"

#' @title Dataset: NCI-60 cancer cell line metrics (k=10 partition)
#' @description
#' Bioinformatics metrics computed on the NCI-60 cancer cell line panel,
#' pre-partitioned into 10 clusters. Used for benchmarking clustering methods.
#' @usage data("nci60_k10")
"nci60_k10"

#' @title Dataset: Golub leukemia gene expression metrics
#' @description
#' Bioinformatics metrics derived from the Golub et al. leukemia gene
#' expression dataset. Used for benchmarking ATSC and clustering methods.
#' @usage data("golub")
"golub"

#' @title Dataset: Breast cancer gene expression metrics
#' @description
#' Bioinformatics metrics derived from a breast cancer gene expression
#' dataset. Used for benchmarking ATSC and PCA-based clustering pipelines.
#' @usage data("breastCancer")
"breastCancer"
