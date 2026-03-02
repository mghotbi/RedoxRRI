#' @title Simulate a Holobiont Redox Dataset with Spatio-Temporal Structure
#'
#' @description
#' Generate a synthetic plant–soil–microbiome dataset with (i) a latent redox
#' state evolving over time, (ii) a mid-season disturbance pulse, (iii) partial
#' cross-domain decoupling, (iv) configurable microbial sparsity (zero inflation),
#' and (v) missing-not-at-random (MNAR) Eh dropout under strongly reduced states.
#'
#' The simulator is designed for benchmarking, validation, and demonstration of
#' RedoxRRI workflows. It is not intended to represent any real ecosystem.
#'
#' @param n_plot Integer >= 1. Number of plots.
#' @param n_depth Integer >= 1. Number of soil depth strata.
#' @param n_plant Integer >= 1. Number of plants per plot–depth.
#' @param n_time Integer >= 2. Number of time points per plant.
#' @param p_micro Integer >= 1. Number of microbial features.
#' @param seed Optional integer seed. Use \code{NULL} for stochastic runs.
#' @param include_graph Logical. If \code{TRUE} and \pkg{igraph} is available,
#'   return a simulated network as \code{graph}.
#' @param depth_labels Optional character vector of length \code{n_depth}.
#'   If \code{NULL}, depths are labeled \code{"D1", "D2", ...}.
#'
#' @param disturbance_strength Numeric in \code{[0, 1]}. Magnitude of disturbance pulse.
#' @param disturbance_center Optional numeric. If \code{NULL}, centered at median time.
#' @param disturbance_width Numeric > 0. Pulse width as fraction of season length.
#'
#' @param seasonal_amp Numeric >= 0. Amplitude of seasonal forcing.
#' @param seasonal_phase Numeric. Phase shift for seasonal sine wave (radians).
#'
#' @param stochastic_reassembly Logical. If \code{TRUE}, microbial counts include a
#'   stochastic reassembly component (mixture of intensities).
#' @param decoupling Numeric in \code{[0, 1]}. Cross-domain decoupling (0=tightly coupled).
#'
#' @param zero_inflation Numeric in \code{[0, 1]}. Probability a microbial entry is set to 0.
#'
#' @param MNAR_strength Numeric in \code{[0, 1]}. Strength of MNAR Eh dropout in reduced states.
#' @param Eh_dropout_threshold Numeric in \code{[0, 1]}. Latent threshold below which dropout can occur.
#'
#' @param micro_mean Numeric > 0. Baseline mean intensity for microbial counts.
#' @param micro_slope Numeric >= 0. Coupling strength between latent redox state and microbial intensity.
#' @param micro_lambda_min Numeric > 0. Lower bound for microbial Poisson intensity.
#' @param micro_lambda_max Numeric > 0. Upper bound for microbial Poisson intensity.
#'
#' @return A named list with:
#' \itemize{
#'   \item \code{id}: Data frame of experimental design (plot, depth, plant_id, time).
#'   \item \code{ROS_flux}: Plant physiological indicators (data frame).
#'   \item \code{Eh_stability}: Soil redox chemistry variables (data frame; Eh may include NA).
#'   \item \code{micro_data}: Microbial feature matrix (data frame; non-negative integers).
#'   \item \code{latent_truth}: Underlying latent redox state in \code{[0, 1]}.
#'   \item \code{graph}: Optional \code{igraph} object if requested and available.
#' }
#'
#' @importFrom stats rnorm rpois runif plogis median
#' @export
simulate_redox_holobiont <- function(
    n_plot = 4,
    n_depth = 2,
    n_plant = 5,
    n_time = 15,
    p_micro = 60,
    seed = 123,
    include_graph = FALSE,
    depth_labels = NULL,
    disturbance_strength = 0.4,
    disturbance_center = NULL,
    disturbance_width = 0.10,
    seasonal_amp = 0.08,
    seasonal_phase = 0,
    stochastic_reassembly = TRUE,
    decoupling = 0.3,
    zero_inflation = 0.4,
    MNAR_strength = 0.4,
    Eh_dropout_threshold = 0.25,
    micro_mean = 8,
    micro_slope = 3,
    micro_lambda_min = 1e-8,
    micro_lambda_max = 1e6
) {
  # ---------------------------
  # helpers
  # ---------------------------
  clamp01 <- function(x) pmin(pmax(x, 0), 1)
  
  .assert_int1 <- function(x, nm, min_val) {
    if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x != as.integer(x) || x < min_val) {
      stop(nm, " must be a single integer >= ", min_val, ".", call. = FALSE)
    }
    invisible(TRUE)
  }
  
  .assert_num1 <- function(x, nm, cond_msg = NULL) {
    if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
      stop(nm, " must be a single finite numeric value.", call. = FALSE)
    }
    if (!is.null(cond_msg)) {
      ok <- tryCatch(eval(parse(text = cond_msg$expr)), error = function(e) FALSE)
      if (!isTRUE(ok)) stop(cond_msg$msg, call. = FALSE)
    }
    invisible(TRUE)
  }
  
  # ---------------------------
  # validation
  # ---------------------------
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
      stop("seed must be NULL or a single finite numeric value.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }
  
  .assert_int1(n_plot,  "n_plot",  1)
  .assert_int1(n_depth, "n_depth", 1)
  .assert_int1(n_plant, "n_plant", 1)
  .assert_int1(n_time,  "n_time",  2)
  .assert_int1(p_micro, "p_micro", 1)
  
  disturbance_strength   <- clamp01(disturbance_strength)
  decoupling             <- clamp01(decoupling)
  zero_inflation         <- clamp01(zero_inflation)
  MNAR_strength          <- clamp01(MNAR_strength)
  Eh_dropout_threshold   <- clamp01(Eh_dropout_threshold)
  
  .assert_num1(disturbance_width, "disturbance_width",
               list(expr = "disturbance_width > 0",
                    msg  = "disturbance_width must be > 0 (fraction of season length)."))
  .assert_num1(seasonal_amp, "seasonal_amp",
               list(expr = "seasonal_amp >= 0",
                    msg  = "seasonal_amp must be >= 0."))
  .assert_num1(seasonal_phase, "seasonal_phase")
  .assert_num1(micro_mean, "micro_mean",
               list(expr = "micro_mean > 0",
                    msg  = "micro_mean must be > 0."))
  .assert_num1(micro_slope, "micro_slope",
               list(expr = "micro_slope >= 0",
                    msg  = "micro_slope must be >= 0."))
  .assert_num1(micro_lambda_min, "micro_lambda_min",
               list(expr = "micro_lambda_min > 0",
                    msg  = "micro_lambda_min must be > 0."))
  .assert_num1(micro_lambda_max, "micro_lambda_max",
               list(expr = "micro_lambda_max > 0",
                    msg  = "micro_lambda_max must be > 0."))
  if (micro_lambda_max <= micro_lambda_min) {
    stop("micro_lambda_max must be > micro_lambda_min.", call. = FALSE)
  }
  
  if (!is.null(depth_labels)) {
    if (!is.character(depth_labels) || length(depth_labels) != n_depth) {
      stop("depth_labels must be NULL or a character vector of length n_depth.", call. = FALSE)
    }
    if (any(!nzchar(depth_labels)) || anyDuplicated(depth_labels)) {
      stop("depth_labels must be non-empty and unique.", call. = FALSE)
    }
  }
  
  # ---------------------------
  # ID structure
  # ---------------------------
  depth_levels <- if (is.null(depth_labels)) paste0("D", seq_len(n_depth)) else depth_labels
  
  id <- expand.grid(
    plot = paste0("P", seq_len(n_plot)),
    depth = depth_levels,
    plant_id = paste0("Plant", seq_len(n_plant)),
    time = seq_len(n_time),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  n <- nrow(id)
  plot_eff <- as.integer(factor(id$plot))
  
  # Depth numeric effect scaled to \code{[0, 1]} (D1 shallow-ish -> 0; deepest -> 1)
  depth_idx <- as.integer(factor(id$depth, levels = depth_levels))
  depth_num <- (depth_idx - 1) / max(n_depth - 1, 1)
  
  tmax <- max(id$time)
  t0 <- if (is.null(disturbance_center)) stats::median(id$time) else as.numeric(disturbance_center)
  if (!is.finite(t0)) t0 <- stats::median(id$time)
  
  sigma2 <- (disturbance_width * tmax)^2
  if (!is.finite(sigma2) || sigma2 <= 0) sigma2 <- (0.10 * tmax)^2
  
  # ---------------------------
  # Latent baseline state (0–1)
  # ---------------------------
  # Baseline: higher depth_num => more reduced => lower latent_truth
  latent_base <- 0.70 -
    0.22 * depth_num +
    stats::rnorm(n_plot, 0, 0.06)[plot_eff]
  
  seasonal <- seasonal_amp * sin(2 * pi * id$time / tmax + seasonal_phase)
  disturbance <- disturbance_strength * exp(-((id$time - t0)^2) / (2 * sigma2))
  
  latent_truth <- latent_base + seasonal - disturbance + stats::rnorm(n, 0, 0.05)
  latent_truth <- clamp01(latent_truth)
  
  # ---------------------------
  # Soil block (Eh–pH coupling + MNAR dropout)
  # ---------------------------
  Eh <- 650 * latent_truth + stats::rnorm(n, 0, 35)
  pH <- 6.5 - 0.8 * (1 - latent_truth) + stats::rnorm(n, 0, 0.2)
  
  Fe2_Fe3 <- stats::plogis(-4 * latent_truth + stats::rnorm(n, 0, 0.6))
  Mn2_Mn4 <- stats::plogis(-3 * latent_truth + stats::rnorm(n, 0, 0.6))
  NH4_NO3 <- stats::plogis(-2.5 * latent_truth + stats::rnorm(n, 0, 0.5))
  
  # MNAR: dropout rises as latent goes below threshold; bounded by MNAR_strength
  denom <- max(Eh_dropout_threshold, 1e-8)
  dropout_prob <- ifelse(
    latent_truth < Eh_dropout_threshold,
    MNAR_strength * (Eh_dropout_threshold - latent_truth) / denom,
    0
  )
  dropout_prob <- clamp01(dropout_prob)
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
  internal_regulation <- decoupling * stats::rnorm(n, 0, 0.10)
  ROS_load <- hypoxia_signal^2 + internal_regulation + stats::rnorm(n, 0, 0.05)
  
  ROS_flux <- data.frame(
    SPAD = 40 - 8 * ROS_load + stats::rnorm(n, 0, 2),
    FvFm = 0.83 - 0.35 * ROS_load + stats::rnorm(n, 0, 0.02),
    PhiPSII = 0.45 - 0.40 * ROS_load + stats::rnorm(n, 0, 0.03),
    NPQ = 0.8 + 2.0 * ROS_load + stats::rnorm(n, 0, 0.2)
  )
  
  # ---------------------------
  # Microbial block (Poisson counts; safe intensities; optional reassembly + sparsity)
  # ---------------------------
  # Intensity increases as system becomes more reduced (low latent_truth)
  # Use log-intensity to avoid overflow, then clamp intensity into [min, max].
  log_lambda <- log(micro_mean) + micro_slope * (1 - latent_truth)
  
  # Prevent pathological exp() overflow/underflow in extreme parameterizations
  log_lambda <- pmin(pmax(log_lambda, log(micro_lambda_min)), log(micro_lambda_max))
  lambda_det <- exp(log_lambda)  # length n
  
  # If stochastic reassembly: mix deterministic intensity with a baseline random intensity
  if (isTRUE(stochastic_reassembly)) {
    # baseline intensity around micro_mean (same scale), independent of latent_truth
    lambda_rand <- rep(micro_mean, n)
    lambda_mix <- (1 - decoupling) * lambda_det + decoupling * lambda_rand
  } else {
    lambda_mix <- lambda_det
  }
  
  lambda_mix <- pmin(pmax(lambda_mix, micro_lambda_min), micro_lambda_max)
  if (any(!is.finite(lambda_mix))) {
    stop("Internal error: microbial intensity became non-finite. Please report with your parameter set.", call. = FALSE)
  }
  
  # Generate counts with per-row intensity replicated across features (fast + stable)
  micro_mat <- matrix(
    stats::rpois(n * p_micro, lambda = rep(lambda_mix, times = p_micro)),
    nrow = n,
    ncol = p_micro,
    byrow = FALSE
  )
  
  # Zero inflation (structural sparsity)
  if (zero_inflation > 0) {
    zero_mask <- matrix(stats::runif(n * p_micro) < zero_inflation, nrow = n, ncol = p_micro)
    micro_mat[zero_mask] <- 0L
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
