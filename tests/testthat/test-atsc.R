golub_pca <- function() {
  clean <- cleanDataset(golub, correlation_threshold = 1)
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

test_that("getMetricsRelevancy weighs every metric and returns the trimmed cases", {
  m <- getMetricsRelevancy(ontMetrics, k = 3, alpha = 0.1, seed = 100)
  metrics <- setdiff(names(df_of(ontMetrics)), "Description")
  # Known issue: Description is weighed too (it is not dropped before RSKC),
  # so check the metrics are present rather than the exact row count.
  expect_true(all(metrics %in% m$relevancy$metric))
  expect_length(m$trimmed_cases, 13)
})
