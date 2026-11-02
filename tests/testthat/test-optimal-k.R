test_that("getOptimalKValue picks the reference k on rnaMetrics", {
  s <- stabilityRange(rnaMetrics, k.range = c(3, 4), bs = 20, seed = 100)
  q <- qualityRange(rnaMetrics, k.range = c(3, 4), seed = 100)
  k <- getOptimalKValue(s, q, k.range = c(3, 4))
  expect_equal(as.numeric(k$Global_optimal_k), c(3, 3))  # returned as character
  expect_equal(k$Stability_max_k_stab, c(0.890138888888889, 1), tolerance = 1e-8)
})

test_that("metricsCorrelations returns the reference Pearson matrix", {
  r <- SummarizedExperiment::assay(metricsCorrelations(rnaMetrics, getImages = FALSE), 1)
  expect_equal(unname(diag(as.matrix(r))), c(1, 1))
  expect_equal(as.matrix(r)["RIN", "DegFact"], -0.974468501362178, tolerance = 1e-8)
})

test_that("plotting functions run without error", {
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plotMetricsMinMax(rnaMetrics))
  expect_no_error(plotMetricsBoxplot(rnaMetrics))
  expect_no_error(plotMetricsCluster(ontMetrics, scale = TRUE))
  expect_no_error(plotMetricsViolin(rnaMetrics))
  expect_no_error(plotMetricsClusterComparison(rnaMetrics, k.vector1 = 3, k.vector2 = c(2, 5)))
})
