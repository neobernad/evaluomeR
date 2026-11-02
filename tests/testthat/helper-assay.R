# Assay of a SummarizedExperiment as a data.frame, one row per metric.
df_of <- function(se) as.data.frame(SummarizedExperiment::assay(se))

# Numeric value of column `col` for metric `metric`.
val <- function(se, metric, col) {
  d <- df_of(se)
  as.numeric(d[d$Metric == metric, col])
}
