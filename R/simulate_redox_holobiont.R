#' Simulate a holobiont redox dataset with spatio-temporal structure
#'
#' @description
#' Generates a synthetic plant–soil–microbe dataset with realistic
#' redox dynamics, disturbance forcing, partial compartment decoupling,
#' and optional stochastic microbial reassembly.
#'
#' Backward compatible with earlier simulator versions.
#'
#' @param n_plot Integer >= 1. Number of plots.
#' @param n_depth Integer (1 or 2). Number of soil depths.
#' @param n_plant Integer >= 1. Plants per plot–depth.
#' @param n_time Integer >= 2. Time points per plant.
#' @param p_micro Integer >= 1. Number of microbial features.
#' @param seed Optional integer seed. NULL for stochastic runs.
#' @param include_graph Logical. If TRUE and igraph is available,
#'   returns a system-level network.
#'
#' @param disturbance_strength Numeric in [0,1].
#'   Controls magnitude of mid-season disturbance pulse.
#' @param stochastic_reassembly Logical.
#'   If TRUE, microbial community includes stochastic reassembly component.
#' @param decoupling Numeric in [0,1].
#'   Degree of cross-domain decoupling (0 = tightly coupled, 1 = decoupled).
#'
#' @return A named list with:
#' \describe{
#'   \item{id}{Sample identifiers.}
#'   \item{ROS_flux}{Plant physiological traits.}
#'   \item{Eh_stability}{Soil redox chemistry variables.}
#'   \item{micro_data}{Zero-inflated microbial abundance matrix.}
#'   \item{latent_truth}{Underlying latent redox state (scaled 0–1).}
#'   \item{graph}{Optional igraph network.}
#' }
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
    include_graph = FALSE,
    disturbance_strength = 0.4,
    stochastic_reassembly = TRUE,
    decoupling = 0.3
) {
  
  # ---------------------------
  # Validation
  # ---------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  if (n_plot < 1 || n_depth < 1 || n_plant < 1 || n_time < 2) {
    stop("Invalid dimensions supplied.")
  }
  
  disturbance_strength <- max(0, min(1, disturbance_strength))
  decoupling <- max(0, min(1, decoupling))
  
  # ---------------------------
  # ID structure
  # ---------------------------
  
  id <- expand.grid(
    plot = paste0("P", seq_len(n_plot)),
    depth = c("shallow", "deep")[seq_len(n_depth)],
    plant_id = paste0("Plant", seq_len(n_plant)),
    time = seq_len(n_time),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  n <- nrow(id)
  
  depth_num <- ifelse(id$depth == "deep", 1, 0)
  plot_eff <- as.numeric(factor(id$plot))
  
  # ---------------------------
  # Latent baseline state
  # ---------------------------
  
  latent_base <- 0.65 -
    0.18 * depth_num +
    stats::rnorm(n_plot, 0, 0.06)[plot_eff]
  
  seasonal <- 0.08 *
    sin(2 * pi * id$time / max(id$time))
  
  disturbance <- disturbance_strength *
    exp(-((id$time - median(id$time))^2) /
          (0.1 * max(id$time)^2))
  
  latent_truth <- latent_base +
    seasonal -
    disturbance +
    stats::rnorm(n, 0, 0.05)
  
  latent_truth <- pmin(pmax(latent_truth, 0), 1)
  
  # ---------------------------
  # Soil block (pH–Eh coupling)
  # ---------------------------
  
  Eh <- 650 * latent_truth +
    stats::rnorm(n, 0, 35)
  
  pH <- 6.5 -
    0.8 * (1 - latent_truth) +
    stats::rnorm(n, 0, 0.2)
  
  Fe2_Fe3 <- stats::plogis(-4 * latent_truth +
                             stats::rnorm(n, 0, 0.6))
  Mn2_Mn4 <- stats::plogis(-3 * latent_truth +
                             stats::rnorm(n, 0, 0.6))
  NH4_NO3 <- stats::plogis(-2.5 * latent_truth +
                             stats::rnorm(n, 0, 0.5))
  
  Eh[latent_truth < 0.25 &
       stats::runif(n) < 0.4] <- NA_real_
  
  Eh_stability <- data.frame(
    Eh = Eh,
    pH = pH,
    Fe2.Fe3 = Fe2_Fe3,
    Mn2.Mn4 = Mn2_Mn4,
    NH4.NO3 = NH4_NO3
  )
  
  # ---------------------------
  # Plant block (partially decoupled)
  # ---------------------------
  
  hypoxia_signal <- (1 - latent_truth)
  
  internal_regulation <- decoupling *
    stats::rnorm(n, 0, 0.1)
  
  ROS_load <- hypoxia_signal^2 +
    internal_regulation +
    stats::rnorm(n, 0, 0.05)
  
  ROS_flux <- data.frame(
    SPAD = 40 - 8 * ROS_load +
      stats::rnorm(n, 0, 2),
    FvFm = 0.83 - 0.35 * ROS_load +
      stats::rnorm(n, 0, 0.02),
    PhiPSII = 0.45 - 0.4 * ROS_load +
      stats::rnorm(n, 0, 0.03),
    NPQ = 0.8 + 2.0 * ROS_load +
      stats::rnorm(n, 0, 0.2)
  )
  
  # ---------------------------
  # Microbial block
  # ---------------------------
  
  lambda <- exp(1.5 - 2 * latent_truth)
  
  micro_mat <- matrix(
    stats::rpois(n * p_micro,
                 lambda = rep(lambda, p_micro)),
    nrow = n
  )
  
  if (isTRUE(stochastic_reassembly)) {
    random_noise <- matrix(
      stats::rpois(n * p_micro, lambda = 1),
      nrow = n
    )
    micro_mat <- (1 - decoupling) * micro_mat +
      decoupling * random_noise
  }
  
  zero_mask <- matrix(
    stats::runif(n * p_micro) < 0.4,
    nrow = n
  )
  micro_mat[zero_mask] <- 0
  
  micro_data <- as.data.frame(micro_mat)
  colnames(micro_data) <- paste0("ASV", seq_len(p_micro))
  
  # ---------------------------
  # Optional network
  # ---------------------------
  
  graph <- NULL
  
  if (isTRUE(include_graph) &&
      requireNamespace("igraph", quietly = TRUE)) {
    
    if (stochastic_reassembly) {
      graph <- igraph::sample_pa(40, m = 2)
    } else {
      graph <- igraph::sample_smallworld(1, 40, 5, 0.1)
    }
    
    graph <- igraph::simplify(graph)
  }
  
  # ---------------------------
  # Return
  # ---------------------------
  
  list(
    id = id,
    ROS_flux = ROS_flux,
    Eh_stability = Eh_stability,
    micro_data = micro_data,
    latent_truth = latent_truth,
    graph = graph
  )
}