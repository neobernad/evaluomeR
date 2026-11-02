# Parallel runs must give the serial result, including when extra arguments
# passed through ... (here kmeansCBI's `runs`) are variables defined in the
# user's session. Those live in the global environment, which is not copied
# to the workers, so the calls are evaluated there on purpose.

run_global <- function(expr) {
  assign("pRuns", 2, envir = globalenv())
  on.exit(rm(pRuns, envir = globalenv()))
  eval(substitute(expr), globalenv())
}

test_that("parallel stability with a global ... argument equals serial", {
  serial <- run_global(stability(rnaMetrics, k = 2, bs = 5, seed = 100, runs = pRuns, numCores = 1))
  par <- run_global(stability(rnaMetrics, k = 2, bs = 5, seed = 100, runs = pRuns, numCores = 2))
  expect_false(anyNA(df_of(par)$Mean_stability_k_2))
  expect_equal(df_of(par), df_of(serial))
})

test_that("parallel quality with a global ... argument equals serial", {
  serial <- run_global(qualityRange(rnaMetrics, k.range = c(2, 3), seed = 100, runs = pRuns, numCores = 1))
  par <- run_global(qualityRange(rnaMetrics, k.range = c(2, 3), seed = 100, runs = pRuns, numCores = 2))
  expect_equal(df_of(par$k_3), df_of(serial$k_3))
})

test_that("one task runs serially even when numCores > 1", {
  expect_no_message(stability(rnaMetrics, k = 2, bs = 5, seed = 100, all_metrics = TRUE,
                              numCores = 2),
                    message = "Number of cores")
})
