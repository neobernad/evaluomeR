test_that("cleanDataset drops non-numeric columns and returns a correlation matrix", {
  r <- cleanDataset(golub, correlation_threshold = 1)
  expect_equal(dim(r$dataset), c(72, 101))
  expect_false("Class" %in% names(r$dataset))
  expect_equal(nrow(r$R), ncol(r$R))
})

test_that("PCASuitability finds ontMetricsOBO unsuitable for PCA", {
  r <- PCASuitability(cleanDataset(ontMetricsOBO, correlation_threshold = 1)$R, sig_level = 0.05)
  expect_false(r$pca_suitable)
})

test_that("performPCA keeps every sample", {
  clean <- cleanDataset(golub, correlation_threshold = 1)$dataset
  set.seed(100) # performPCA's nFactors bootstrap takes no seed
  p <- performPCA(clean)
  expect_equal(dim(p$dataset_ncp), c(72, 20))
})
