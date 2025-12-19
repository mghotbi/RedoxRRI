testthat::test_that("RRI_pipeline returns expected structure and ranges", {
  data("redoxrri_example", package = "RedoxRRI")

  out <- RedoxRRI::RRI_pipeline(
    redoxrri_example$ROS_flux,
    redoxrri_example$Eh_stability,
    micro_data  = redoxrri_example$micro_data,
    graph       = redoxrri_example$graph,
    alpha_micro = 0.6,
    method_phys = "pca",
    method_soil = "pca",
    method_micro = "pca"
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_true(all(c("Physio", "Soil", "Micro", "RRI", "RRI_index") %in% names(out)))

  # compositional rows should sum to 1 (allow tiny numeric tolerance)
  s <- rowSums(out[, c("Physio", "Soil", "Micro")])
  testthat::expect_true(all(abs(s - 1) < 1e-8))

  # RRI should be in [0, 1]
  testthat::expect_true(all(out$RRI >= 0 & out$RRI <= 1, na.rm = TRUE))

  # attribute exists
  testthat::expect_true(is.numeric(attr(out, "RRI_index")))
})
