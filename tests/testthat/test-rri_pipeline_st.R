library(testthat)

test_that("rri_pipeline_st returns a valid RRI object (snapshot)", {
  ROS_flux <- data.frame(
    FvFm = c(0.8, 0.7, 0.75, 0.72),
    CAT  = c(10, 12, 11, 9)
  )

  Eh_stability <- data.frame(
    Eh = c(350, 320, 340, 310),
    pH = c(6.5, 6.7, 6.6, 6.8)
  )

  micro_data <- data.frame(
    geneA = c(1, 2, 1, 3),
    geneB = c(5, 4, 6, 3)
  )

  id <- data.frame(
    plot  = c("A", "A", "B", "B"),
    depth = c("shallow", "deep", "shallow", "deep"),
    time  = 1:4
  )

  res <- rri_pipeline_st(
    ROS_flux = ROS_flux,
    Eh_stability = Eh_stability,
    micro_data = micro_data,
    id = id,
    group_cols = c("plot", "depth"),
    scale_by = c("plot", "depth"),
    direction_phys = "auto",
    direction_anchor_phys = "FvFm",
    direction_soil = "auto",
    direction_anchor_soil = "Eh"
  )

  expect_s3_class(res, "RRI")
  expect_named(res, c("row_scores", "row_scores_comp", "dyn_scores", "meta"))

  expect_s3_class(res$row_scores, "data.frame")
  expect_equal(nrow(res$row_scores), nrow(ROS_flux))

  expect_true(all(c("Physio", "Soil", "Micro", "RRI") %in% names(res$row_scores)))
  expect_true(all(res$row_scores$RRI >= 0 & res$row_scores$RRI <= 1, na.rm = TRUE))

  expect_null(res$dyn_scores)

  expect_true(is.list(res$meta))
  expect_true(is.numeric(res$meta$rri_index))
})


test_that("rri_pipeline_st enforces compositional scaling", {
  ROS_flux <- data.frame(a = 1:3, b = 2:4)
  Eh_stability <- data.frame(c = 3:5, d = 6:8)
  micro_data <- data.frame(e = 1:3)

  res <- rri_pipeline_st(
    ROS_flux = ROS_flux,
    Eh_stability = Eh_stability,
    micro_data = micro_data
  )

  comp <- res$row_scores_comp[, c("Physio", "Soil", "Micro")]
  rs <- rowSums(comp)

  expect_true(all(abs(rs[is.finite(rs)] - 1) < 1e-8))
})


test_that("rri_pipeline_st errors when microbial input is missing", {
  ROS_flux <- data.frame(a = 1:3)
  Eh_stability <- data.frame(b = 2:4)

  expect_error(
    rri_pipeline_st(
      ROS_flux = ROS_flux,
      Eh_stability = Eh_stability
    ),
    "Provide microbial input"
  )
})
