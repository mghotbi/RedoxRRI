#' @title Sensitivity Analysis for RedoxRRI
#'
#' @description
#' Evaluates robustness of RRI rankings under alternative domain
#' aggregation weights using already computed domain scores.
#'
#' @param res An object returned by \code{rri_pipeline_st()}.
#' @param weight_grid Numeric vector of Physio weights in (0,1).
#'
#' @return A data frame containing tested weights and Spearman rank
#'   correlations relative to the baseline RRI.
#'
#' @details
#' This function perturbs only the aggregation step and does not
#' recompute latent domain scores. It therefore evaluates stability
#' of domain integration rather than reducer sensitivity.
#' @examples
#' # ---- Simulate a small holobiont dataset ----
#' sim <- simulate_redox_holobiont(
#'   n_plot = 2,
#'   n_depth = 4,
#'   n_plant = 2,
#'   n_time = 8,
#'   p_micro = 20,
#'   seed = 1
#' )
#'
#' # ---- Compute RedoxRRI ----
#' res <- rri_pipeline_st(
#'   ROS_flux = sim$ROS_flux,
#'   Eh_stability = sim$Eh_stability,
#'   micro_data = sim$micro_data,
#'   id = sim$id,
#'   reducer = "per_domain",
#'   scaling = "pnorm"
#' )
#'
#' # ---- Evaluate aggregation sensitivity ----
#' sens <- rri_sensitivity(
#'   res,
#'   weight_grid = seq(0.3, 0.5, by = 0.1)
#' )
#'
#' 
#' @export
rri_sensitivity <- function(res,
                            weight_grid = seq(0.2, 0.6, by = 0.1)) {
  
  required <- c("Physio", "Soil", "Micro", "RRI")
  
  if (!all(required %in% names(res$row_scores))) {
    stop("res$row_scores must contain Physio, Soil, Micro, and RRI.",
         call. = FALSE)
  }
  
  if (!is.numeric(weight_grid) || any(!is.finite(weight_grid))) {
    stop("weight_grid must be a numeric vector of finite values.",
         call. = FALSE)
  }
  
  if (any(weight_grid <= 0 | weight_grid >= 1)) {
    stop("All weights must lie strictly between 0 and 1.",
         call. = FALSE)
  }
  
  dom <- res$row_scores[, c("Physio", "Soil", "Micro")]
  baseline_rri <- res$row_scores$RRI
  
  results <- lapply(weight_grid, function(w_phys) {
    
    remaining <- 1 - w_phys
    w_soil  <- remaining / 2
    w_micro <- remaining / 2
    
    alt_rri <- w_phys  * dom$Physio +
      w_soil  * dom$Soil +
      w_micro * dom$Micro
    
    cor_val <- stats::cor(
      baseline_rri,
      alt_rri,
      method = "spearman",
      use = "complete.obs"
    )
    
    data.frame(
      weight_physio = w_phys,
      weight_soil   = w_soil,
      weight_micro  = w_micro,
      spearman_rank_correlation = cor_val
    )
  })
  
  do.call(rbind, results)
}