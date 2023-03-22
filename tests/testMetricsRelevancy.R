library(evaluomeR)

data("ontMetrics")
metricsRelevancy = getMetricsRelevancy(ontMetrics, k=3, alpha=0.1, seed=100)
# RSKC output object
metricsRelevancy$rskc
# Trimmed cases from input (row indexes)
metricsRelevancy$trimmed_cases
# Metrics relevancy table
metricsRelevancy$relevancy
