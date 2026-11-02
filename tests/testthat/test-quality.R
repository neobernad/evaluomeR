test_that("quality matches reference silhouette widths on rnaMetrics", {
  q <- quality(rnaMetrics, k = 3, seed = 13606)
  expect_equal(val(q, "RIN", "Avg_Silhouette_Width"), 0.627829396038413, tolerance = 1e-8)
  expect_equal(val(q, "DegFact", "Avg_Silhouette_Width"), 0.737191191352892, tolerance = 1e-8)
  expect_equal(val(q, "RIN", "Cluster_1_Size") + val(q, "RIN", "Cluster_2_Size") +
                 val(q, "RIN", "Cluster_3_Size"), nrow(df_of(rnaMetrics)))
})

test_that("qualityRange matches reference values and getDataQualityRange picks a k", {
  qr <- qualityRange(rnaMetrics, k.range = c(2, 4), seed = 20)
  expect_named(qr, c("k_2", "k_3", "k_4"))
  q2 <- getDataQualityRange(qr, 2)
  q4 <- getDataQualityRange(qr, 4)
  expect_equal(val(q2, "RIN", "Avg_Silhouette_Width"), 0.608402004052639, tolerance = 1e-8)
  expect_equal(val(q2, "DegFact", "Avg_Silhouette_Width"), 0.666587617027136, tolerance = 1e-8)
  expect_equal(val(q4, "RIN", "Avg_Silhouette_Width"), 0.463905611516569, tolerance = 1e-8)
  expect_equal(val(q4, "DegFact", "Avg_Silhouette_Width"), 0.634170498361632, tolerance = 1e-8)
})

test_that("qualitySet returns one entry per k", {
  qs <- qualitySet(rnaMetrics, k.set = c(2, 4), seed = 13606)
  expect_named(qs, c("k_2", "k_4"))
})

test_that("Calinski-Harabasz index matches reference values", {
  ch <- quality(ontMetrics, k = 3, quality_index = "ch", seed = 100)
  expect_equal(nrow(ch), 19)
  expect_equal(mean(ch$CH), 796.444197458028, tolerance = 1e-8)
  expect_equal(ch$CH[ch$Metric == "ANOnto"], 319.846916590098, tolerance = 1e-8)

  chr <- qualityRange(ontMetrics, k.range = c(2, 3), quality_index = "ch", seed = 100)
  expect_length(chr, 2)
})

test_that("RSKC quality runs with all metrics merged", {
  q <- qualityRange(ontMetrics, k.range = c(3, 3), seed = 13007, all_metrics = TRUE,
                    cbi = "rskc", L1 = 2, alpha = 0)
  expect_equal(val(q$k_3, "all_metrics", "Cluster_1_SilScore"), 0.35653231531208, tolerance = 1e-8)
})
