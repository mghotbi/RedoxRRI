testthat::test_that("plot_RRI_ternary errors if required columns are missing", {
  testthat::skip_if_not_installed("ggtern")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("viridis")
  testthat::skip_if_not_installed("igraph")

  bad <- data.frame(a = 1:3)
  testthat::expect_error(plot_RRI_ternary(bad), "missing required columns")
})

testthat::test_that("plot_RRI_ternary returns a ggplot when suggested packages exist", {
  testthat::skip_if_not_installed("ggtern")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("viridis")

  data("redoxrri_example", package = "RedoxRRI")

  out <- RedoxRRI::RRI_pipeline(
    redoxrri_example$ROS_flux,
    redoxrri_example$Eh_stability,
    micro_data = redoxrri_example$micro_data,
    graph      = redoxrri_example$graph
  )

  p <- RedoxRRI::plot_RRI_ternary(out)
  testthat::expect_s3_class(p, "ggplot")
})
