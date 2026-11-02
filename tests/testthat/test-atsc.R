# performPCA runs an nFactors bootstrap that takes no seed and crashes for a
# few RNG states on golub (e.g. set.seed(17)), so seed it here.
golub_pca <- function() {
  clean <- cleanDataset(golub, correlation_threshold = 1)
  set.seed(100)
  performPCA(clean$dataset)$dataset_ncp
}

test_that("ATSC finds k = 3 on golub and trims the reference rows", {
  r <- ATSC(golub_pca(), alpha = 0.1, k.range = c(3, 10), cbi = "kmeans")
  expect_equal(r$optimalK, 3)
  expect_equal(r$trimmedRows, c(2, 17, 18, 25, 30, 31, 32, 35, 38, 69))
  expect_equal(nrow(r$trimmmedDataset), 72 - 10)
})

test_that("ATSC trims the gold standard together with the rows", {
  gs <- as.numeric(factor(golub$Class, levels = c("B", "T", "M")))
  r <- ATSC(golub_pca(), alpha = 0.1, k.range = c(3, 10), cbi = "kmeans", gold_standard = gs)
  expect_equal(r$optimalK, 3)
  expect_length(r$gold_standard_trimmed, 62)
})

test_that("getRSKCL1Boundry is reproducible whatever the RNG state", {
  set.seed(1)
  a <- getRSKCL1Boundry(ontMetrics, k = 3, seed = 100)
  set.seed(2)
  b <- getRSKCL1Boundry(ontMetrics, k = 3, seed = 100)
  expect_identical(a, b)
  expect_equal(a, 2)
})

test_that("getMetricsRelevancy weighs only the metrics, reproducibly", {
  m <- getMetricsRelevancy(ontMetrics, k = 3, alpha = 0.1, seed = 100)
  metrics <- setdiff(names(df_of(ontMetrics)), "Description")
  expect_setequal(m$relevancy$metric, metrics)
  expect_equal(m$relevancy$metric[1], "WMCOnto")
  expect_equal(m$relevancy$weight[1], 0.99267148870117661, tolerance = 1e-8)
  expect_equal(m$trimmed_cases, c(2, 3, 5, 13, 21, 41, 59, 67, 71, 73))
  expect_identical(getMetricsRelevancy(ontMetrics, k = 3, alpha = 0.1, seed = 100)$relevancy, m$relevancy)
})
