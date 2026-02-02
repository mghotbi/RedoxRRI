#' Simulate a holobiont redox dataset with spatio-temporal structure
#'
#' @description
#' Generates a synthetic plant–soil–microbe dataset with realistic properties
#' commonly encountered in redox ecology and holobiont studies.
#'
#' @param n_plot Integer >= 1. Number of plots.
#' @param n_depth Integer equal to 1 or 2. Number of soil depths (shallow/deep).
#' @param n_plant Integer >= 1. Number of plants per plot–depth.
#' @param n_time Integer >= 2. Number of time points per plant.
#' @param p_micro Integer >= 1. Number of microbial features (e.g. ASVs).
#' @param seed Optional integer seed; use \code{NULL} for stochastic simulation.
#' @param include_graph Logical. If TRUE and \pkg{igraph} is available, returns a network.
#'
#' @return A named list with components:
#' \describe{
#'   \item{id}{data.frame of sample identifiers: plot, depth, plant_id, time.}
#'   \item{ROS_flux}{data.frame of plant physiological traits.}
#'   \item{Eh_stability}{data.frame of soil redox chemistry variables.}
#'   \item{micro_data}{data.frame of zero-inflated microbial abundances.}
#'   \item{latent_truth}{Numeric vector scaled to \eqn{[0,1]} representing latent redox state.}
#'   \item{graph}{Optional \pkg{igraph} object (system-level network).}
#' }
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_redox_holobiont(seed = 1)
#' str(sim)
#'
#' @importFrom stats rnorm rpois runif plogis
#' @export
simulate_redox_holobiont <- function(
  n_plot = 4,
  n_depth = 2,
  n_plant = 5,
  n_time = 15,
  p_micro = 60,
  seed = 123,
  include_graph = FALSE
) {
  id <- expand.grid(
    plot = paste0("P", seq_len(n_plot)),
    depth = c("shallow", "deep")[seq_len(n_depth)],
    plant_id = paste0("Plant", seq_len(n_plant)),
    time = seq_len(n_time)
  )

  n <- nrow(id)
  depth_num <- ifelse(id$depth == "deep", 1, 0)
  plot_eff <- as.numeric(factor(id$plot))

  latent_truth <- 0.65 -
    0.18 * depth_num +
    stats::rnorm(n_plot, 0, 0.06)[plot_eff] +
    0.08 * sin(2 * pi * id$time / max(id$time)) +
    stats::rnorm(n, 0, 0.05)

  latent_truth <- pmin(pmax(latent_truth, 0), 1)

  Eh <- 650 * latent_truth + stats::rnorm(n, 0, 35)
  Fe2_Fe3 <- stats::plogis(-4 * latent_truth + stats::rnorm(n, 0, 0.6))
  Mn2_Mn4 <- stats::plogis(-3 * latent_truth + stats::rnorm(n, 0, 0.6))
  NH4_NO3 <- stats::plogis(-2.5 * latent_truth + stats::rnorm(n, 0, 0.5))

  Eh[latent_truth < 0.25 & stats::runif(n) < 0.4] <- NA

  Eh_stability <- data.frame(
    Eh = Eh,
    Fe2.Fe3 = Fe2_Fe3,
    Mn2.Mn4 = Mn2_Mn4,
    NH4.NO3 = NH4_NO3
  )

  ROS_load <- (1 - latent_truth)^2 + stats::rnorm(n, 0, 0.05)

  ROS_flux <- data.frame(
    SPAD = 40 - 8 * ROS_load + stats::rnorm(n, 0, 2),
    FvFm = 0.83 - 0.35 * ROS_load + stats::rnorm(n, 0, 0.02),
    PhiPSII = 0.45 - 0.4 * ROS_load + stats::rnorm(n, 0, 0.03),
    NPQ = 0.8 + 2.0 * ROS_load + stats::rnorm(n, 0, 0.2)
  )

  lambda <- exp(1.5 - 2 * latent_truth)

  micro_mat <- matrix(
    stats::rpois(n * p_micro, lambda = rep(lambda, p_micro)),
    nrow = n
  )

  zero_mask <- matrix(stats::runif(n * p_micro) < 0.4, nrow = n)
  micro_mat[zero_mask] <- 0

  micro_data <- as.data.frame(micro_mat)
  colnames(micro_data) <- paste0("ASV", seq_len(p_micro))

  graph <- NULL
  if (isTRUE(include_graph) && requireNamespace("igraph", quietly = TRUE)) {
    graph <- igraph::sample_smallworld(1, 40, 5, 0.1)
    graph <- igraph::simplify(graph)
  }

  list(
    id = id,
    ROS_flux = ROS_flux,
    Eh_stability = Eh_stability,
    micro_data = micro_data,
    latent_truth = latent_truth,
    graph = graph
  )
}
