test_that("cluster-resolution ranking prefers the highest resolution", {
  hits <- c(
    "integratedRPCA_snn_res.0.8",
    "integratedRPCA_snn_res.3.0",
    "integratedRPCA_snn_res.1.5"
  )
  expect_equal(
    uenoyscRNA:::cf_rank_pattern_hits(hits, role = "feature")[[1]],
    "integratedRPCA_snn_res.3.0"
  )
})

test_that("natural level ordering handles numeric cluster labels", {
  expect_equal(
    uenoyscRNA:::cf_natural_levels(c("10", "2", "1")),
    c("1", "2", "10")
  )
})
