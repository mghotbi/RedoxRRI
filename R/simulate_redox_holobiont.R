#' @title Simulate a Holobiont Redox Dataset with Spatio-Temporal Structure
#'
#' @description
#' Generates a synthetic plant–soil–microbe dataset with structured redox
#' dynamics, disturbance forcing, partial cross-domain decoupling, configurable
#' sparsity (zero inflation), and missing-not-at-random (MNAR) Eh dropout.
#' The simulator is intended for benchmarking, methodological validation,
#' and demonstration of RedoxRRI workflows.
#'
#' @param n_plot Integer greater than or equal to 1. Number of plots.
#' @param n_depth Integer (1 or 2). Number of soil depths.
#' @param n_plant Integer greater than or equal to 1. Number of plants per plot–depth.
#' @param n_time Integer greater than or equal to 2. Number of time points per plant.
#' @param p_micro Integer greater than or equal to 1. Number of microbial features.
#' @param seed Optional integer seed. Use \code{NULL} for stochastic runs.
#' @param include_graph Logical. If \code{TRUE} and \pkg{igraph} is available,
#'   returns a simulated system-level network.
#'
#' @param disturbance_strength Numeric between 0 and 1. Magnitude of the
#'   mid-season disturbance pulse.
#' @param disturbance_center Optional numeric. If \code{NULL}, the disturbance
#'   is centred at the median time point.
#' @param disturbance_width Numeric greater than 0. Controls pulse width
#'   as a fraction of the total season length.
#'
#' @param seasonal_amp Numeric greater than or equal to 0. Amplitude of
#'   seasonal forcing.
#' @param seasonal_phase Numeric. Phase shift for the seasonal sine wave (radians).
#'
#' @param stochastic_reassembly Logical. If \code{TRUE}, microbial communities
#'   combine deterministic structure with a stochastic reassembly component.
#' @param decoupling Numeric between 0 and 1 indicating cross-domain
#'   decoupling (0 = tightly coupled, 1 = fully decoupled).
#'
#' @param zero_inflation Numeric between 0 and 1 giving the probability that
#'   a microbial entry is set to zero.
#'
#' @param MNAR_strength Numeric between 0 and 1 controlling the strength of
#'   missing-not-at-random Eh dropout under strongly reduced states.
#' @param Eh_dropout_threshold Numeric between 0 and 1 specifying the latent
#'   redox threshold below which Eh dropout may occur.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{id}: Experimental design structure.
#'   \item \code{ROS_flux}: Simulated plant physiological indicators.
#'   \item \code{Eh_stability}: Soil redox chemistry variables.
#'   \item \code{micro_data}: Microbial feature matrix.
#'   \item \code{latent_truth}: Simulated underlying redox state scaled between 0 and 1.
#'   \item \code{graph}: Optional \code{igraph} object if requested.
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
    disturbance_center = NULL,
    disturbance_width = 0.10,
    seasonal_amp = 0.08,
    seasonal_phase = 0,
    stochastic_reassembly = TRUE,
    decoupling = 0.3,
    zero_inflation = 0.4,
    MNAR_strength = 0.4,
    Eh_dropout_threshold = 0.25
) {
  # ---------------------------
  # Validation
  # ---------------------------
  if (!is.null(seed)) set.seed(seed)
  
  stopifnot(n_plot >= 1, n_depth >= 1, n_plant >= 1, n_time >= 2, p_micro >= 1)
  
  clamp01 <- function(x) max(0, min(1, x))
  disturbance_strength <- clamp01(disturbance_strength)
  decoupling <- clamp01(decoupling)
  zero_inflation <- clamp01(zero_inflation)
  MNAR_strength <- clamp01(MNAR_strength)
  Eh_dropout_threshold <- clamp01(Eh_dropout_threshold)
  
  if (!is.numeric(disturbance_width) || length(disturbance_width) != 1 || disturbance_width <= 0) {
    stop("disturbance_width must be a single numeric > 0 (interpreted as fraction of season length).")
  }
  if (!is.numeric(seasonal_amp) || seasonal_amp < 0) stop("seasonal_amp must be >= 0.")
  
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
  
  tmax <- max(id$time)
  t0 <- if (is.null(disturbance_center)) stats::median(id$time) else as.numeric(disturbance_center)
  sigma2 <- (disturbance_width * tmax)^2
  if (!is.finite(sigma2) || sigma2 <= 0) sigma2 <- (0.10 * tmax)^2
  
  # ---------------------------
  # Latent baseline state (0–1)
  # ---------------------------
  latent_base <- 0.65 -
    0.18 * depth_num +
    stats::rnorm(n_plot, 0, 0.06)[plot_eff]
  
  seasonal <- seasonal_amp * sin(2 * pi * id$time / tmax + seasonal_phase)
  
  disturbance <- disturbance_strength * exp(-((id$time - t0)^2) / (2 * sigma2))
  
  latent_truth <- latent_base + seasonal - disturbance + stats::rnorm(n, 0, 0.05)
  latent_truth <- pmin(pmax(latent_truth, 0), 1)
  
  # ---------------------------
  # Soil block (pH–Eh coupling + MNAR dropout)
  # ---------------------------
  Eh <- 650 * latent_truth + stats::rnorm(n, 0, 35)
  
  pH <- 6.5 - 0.8 * (1 - latent_truth) + stats::rnorm(n, 0, 0.2)
  
  Fe2_Fe3 <- stats::plogis(-4 * latent_truth + stats::rnorm(n, 0, 0.6))
  Mn2_Mn4 <- stats::plogis(-3 * latent_truth + stats::rnorm(n, 0, 0.6))
  NH4_NO3 <- stats::plogis(-2.5 * latent_truth + stats::rnorm(n, 0, 0.5))
  
  # MNAR Eh dropout increases as latent_truth falls below threshold
  # dropout_prob in [0, MNAR_strength], rising linearly as latent decreases
  dropout_prob <- ifelse(
    latent_truth < Eh_dropout_threshold,
    MNAR_strength * (Eh_dropout_threshold - latent_truth) / max(Eh_dropout_threshold, 1e-8),
    0
  )
  Eh[stats::runif(n) < dropout_prob] <- NA_real_
  
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
  
  internal_regulation <- decoupling * stats::rnorm(n, 0, 0.1)
  
  ROS_load <- hypoxia_signal^2 + internal_regulation + stats::rnorm(n, 0, 0.05)
  
  ROS_flux <- data.frame(
    SPAD = 40 - 8 * ROS_load + stats::rnorm(n, 0, 2),
    FvFm = 0.83 - 0.35 * ROS_load + stats::rnorm(n, 0, 0.02),
    PhiPSII = 0.45 - 0.4 * ROS_load + stats::rnorm(n, 0, 0.03),
    NPQ = 0.8 + 2.0 * ROS_load + stats::rnorm(n, 0, 0.2)
  )
  
  # ---------------------------
  # Microbial block (counts + configurable zero inflation)
  # ---------------------------
  lambda <- exp(1.5 - 2 * latent_truth)
  
  micro_mat <- matrix(
    stats::rpois(n * p_micro, lambda = rep(lambda, p_micro)),
    nrow = n
  )
  
  if (isTRUE(stochastic_reassembly)) {
    random_noise <- matrix(stats::rpois(n * p_micro, lambda = 1), nrow = n)
    micro_mat <- (1 - decoupling) * micro_mat + decoupling * random_noise
  }
  
  if (zero_inflation > 0) {
    zero_mask <- matrix(stats::runif(n * p_micro) < zero_inflation, nrow = n)
    micro_mat[zero_mask] <- 0
  }
  
  micro_data <- as.data.frame(micro_mat)
  colnames(micro_data) <- paste0("ASV", seq_len(p_micro))
  
  # ---------------------------
  # Optional network
  # ---------------------------
  graph <- NULL
  if (isTRUE(include_graph) && requireNamespace("igraph", quietly = TRUE)) {
    if (isTRUE(stochastic_reassembly)) {
      graph <- igraph::sample_pa(40, m = 2)
    } else {
      graph <- igraph::sample_smallworld(1, 40, 5, 0.1)
    }
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