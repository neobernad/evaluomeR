# Reference values come from the package with the seeds given. A change in
# them means the bootstrap or clustering changed, not just a refactor.

test_that("stabilityRange matches reference values on rnaMetrics", {
  s <- stabilityRange(rnaMetrics, k.range = c(2, 4), bs = 20, seed = 13606)
  expect_equal(val(s, "RIN", "Mean_stability_k_2"), 0.825833333333333, tolerance = 1e-8)
  expect_equal(val(s, "RIN", "Mean_stability_k_3"), 0.778412698412698, tolerance = 1e-8)
  expect_equal(val(s, "RIN", "Mean_stability_k_4"), 0.69625, tolerance = 1e-8)
  expect_equal(val(s, "DegFact", "Mean_stability_k_2"), 0.955595238095238, tolerance = 1e-8)
  expect_equal(val(s, "DegFact", "Mean_stability_k_3"), 0.977777777777778, tolerance = 1e-8)
  expect_equal(val(s, "DegFact", "Mean_stability_k_4"), 0.820833333333333, tolerance = 1e-8)
})

test_that("stability for one k equals the same k inside stabilityRange", {
  s <- stability(rnaMetrics, k = 4, bs = 20, seed = 13606)
  expect_equal(val(s, "RIN", "Mean_stability_k_4"), 0.69625, tolerance = 1e-8)
  expect_equal(val(s, "DegFact", "Mean_stability_k_4"), 0.820833333333333, tolerance = 1e-8)
})

test_that("stabilitySet returns one stability column per k", {
  s <- stabilitySet(rnaMetrics, k.set = c(2, 4), bs = 20, seed = 13606)
  d <- df_of(s)
  expect_setequal(d$Metric, c("RIN", "DegFact"))
  expect_true(all(c("Mean_stability_k_2", "Mean_stability_k_4") %in% names(d)))
})

test_that("all_metrics stability equals clusterbootWrapper on the same data", {
  x <- as.data.frame(SummarizedExperiment::assay(rnaMetrics))
  s <- stability(x, k = 2, bs = 100, seed = 100, all_metrics = TRUE)
  merged <- val(s$stability_mean, "all_metrics", "Mean_stability_k_2")
  boot <- clusterbootWrapper(data = x[c("RIN", "DegFact")], B = 100, bootmethod = "boot",
                             cbi = "kmeans", krange = 2, seed = 100, gold_standard = NULL)
  expect_equal(mean(boot$bootmean), merged, tolerance = 1e-8)
})

test_that("every supported CBI gives its reference stability", {
  expected <- c(kmeans = 0.807202380952381, clara = 0.897261904761905,
                clara_pam = 0.849047619047619, hclust = 0.813630952380952,
                pamk = 0.897261904761905, pamk_pam = 0.849047619047619)
  for (cbi in names(expected)) {
    s <- stability(rnaMetrics, cbi = cbi, k = 2, bs = 10, seed = 13606)
    expect_equal(val(s, "RIN", "Mean_stability_k_2"), expected[[cbi]], tolerance = 1e-8, label = cbi)
  }
  expect_true(all(names(expected) %in% evaluomeRSupportedCBI()))
})

test_that("a gold standard gives ARI per k on golub", {
  gs <- as.numeric(factor(golub$Class, levels = c("B", "T", "M")))
  clean <- cleanDataset(golub, correlation_threshold = 1)$dataset
  s <- stabilityRange(clean, k.range = c(2, 3), bs = 10, all_metrics = TRUE, seed = 100, gold_standard = gs)
  ari <- standardizeStabilityData(s, k.range = c(2, 3))
  expect_equal(unlist(ari["all_metrics", ]), c(k_2 = 0.693041608535879, k_3 = 0.913043422965347),
               tolerance = 1e-8)
})

test_that("a two-class gold standard gives ARI per k on breastCancer", {
  gs <- as.numeric(factor(breastCancer$diagnosis, levels = c("B", "M")))
  s <- stabilityRange(breastCancer[, -2], k.range = c(2, 3), bs = 10, all_metrics = TRUE,
                      seed = 100, gold_standard = gs)
  ari <- standardizeStabilityData(s, k.range = c(2, 3))
  expect_equal(unlist(ari["all_metrics", ]), c(k_2 = 0.491424536224555, k_3 = 0.541206168253931),
               tolerance = 1e-8)
})
