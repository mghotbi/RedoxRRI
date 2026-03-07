#' @title Example workflow for RedoxRRI
#'
#' @name redoxrri_example
#'
#' @description
#' Demonstrates the complete workflow for computing the
#' Holobiont Redox Resilience Index (RRI) using simulated
#' plant, soil, and microbiome data.
#'
#' @details
#' This example generates a fully synthetic dataset using
#' \code{\link{simulate_redox_holobiont}} and then computes the
#' Redox Resilience Index using \code{\link{rri_pipeline_st}}.
#'
#' The simulation framework represents a simplified
#' soil–plant–microbiome system and includes:
#'
#' \itemize{
#'   \item Disturbance forcing over time
#'   \item Partial plant–soil–microbiome decoupling
#'   \item pH–Eh interaction within the soil redox domain
#'   \item Zero-inflated high-dimensional microbial abundance data
#'   \item Optional stochastic microbial reassembly
#'   \item A microbial association network generated using \pkg{igraph}
#' }
#'
#' The simulated data do \emph{not} represent any real biological
#' system and are provided solely for demonstration purposes.
#'
#' @examples
#' # Generate a synthetic dataset
#' sim <- simulate_redox_holobiont(seed = 1)
#'
#' # Inspect structure
#' str(sim)
#'
#' # Compute the Redox Resilience Index
#' res <- rri_pipeline_st(
#'   ROS_flux = sim$ROS_flux,
#'   Eh_stability = sim$Eh_stability,
#'   micro_data = sim$micro_data,
#'   graph = sim$graph,
#'   id = sim$id,
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
#' @seealso
#' \code{\link{simulate_redox_holobiont}},
#' \code{\link{rri_pipeline_st}},
#' \code{\link{plot_RRI_ternary}}
#'
NULL