#' Example data for RedoxRRI
#'
#' A small synthetic dataset bundled with \pkg{RedoxRRI} to demonstrate
#' the complete workflow for computing the Holobiont Redox Resilience
#' Index (RRI) and visualizing plant–soil–microbe contributions.
#'
#' @details
#' This dataset is \strong{fully synthetic} and generated using
#' \code{\link{simulate_redox_holobiont}}.
#'
#' It is intended exclusively for:
#' \itemize{
#'   \item Function examples
#'   \item Unit tests
#'   \item Vignettes and tutorials
#' }
#'
#' The simulation includes:
#' \itemize{
#'   \item Disturbance forcing over time
#'   \item Partial plant–soil–microbe decoupling
#'   \item pH–Eh interaction in the soil block
#'   \item Zero-inflated high-dimensional microbial data
#'   \item Optional stochastic microbial reassembly
#'   \item A system-level microbial association network (if igraph available)
#' }
#'
#' The data do \emph{not} represent any real biological system.
#'
#' @format A named list with components:
#' \describe{
#'   \item{id}{data.frame; sample identifiers (plot, depth, plant_id, time).}
#'   \item{ROS_flux}{data.frame; plant physiological traits.}
#'   \item{Eh_stability}{data.frame; soil redox chemistry variables.}
#'   \item{micro_data}{data.frame; sparse microbial abundance matrix.}
#'   \item{latent_truth}{numeric vector; underlying simulated redox state scaled to [0,1].}
#'   \item{graph}{Optional \pkg{igraph} object representing a microbial association network.}
#' }
#'
#' @source Generated via \code{simulate_redox_holobiont(seed = 1)}.
#'
#' @examples
#' data("redoxrri_example")
#'
#' str(redoxrri_example)
#'
#' res <- rri_pipeline_st(
#'   ROS_flux = redoxrri_example$ROS_flux,
#'   Eh_stability = redoxrri_example$Eh_stability,
#'   micro_data = redoxrri_example$micro_data,
#'   graph = redoxrri_example$graph,
#'   id = redoxrri_example$id,
#'   group_cols = c("plot", "depth"),
#'   scale_by = c("plot", "depth"),
#'   direction_phys = "auto",
#'   direction_anchor_phys = "FvFm",
#'   direction_soil = "auto",
#'   direction_anchor_soil = "Eh"
#' )
#'
#' head(res$row_scores)
#' head(res$row_scores_comp)
#'
#' @docType data
#' @name redoxrri_example
NULL