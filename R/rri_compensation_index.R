#' @title Cross-Domain Compensation Index
#'
#' @description
#' Computes a scale-invariant cross-domain compensation index
#' based on mean pairwise correlations among Physio, Soil,
#' and Micro domain scores.
#'
#' @param res An object returned by \code{rri_pipeline_st()}.
#'
#' @return A single numeric compensation index in \eqn{[-1, 1]}.
#'   Positive values indicate compensatory dynamics
#'   (negative cross-domain correlations),
#'   while negative values indicate synchrony.
#'
#' @details
#' The compensation index is defined as:
#'
#' \deqn{
#'   \mathrm{Comp} = -\frac{2}{k(k-1)} \sum_{i<j} \mathrm{Cor}(D_i, D_j)
#' }
#'
#' where \eqn{k = 3} and \eqn{D_i} are domain-level scores.
#'
#' Using correlations rather than covariances ensures scale invariance
#' and boundedness across datasets.
#' @examples
#' # ---- Simulate small holobiont dataset ----
#' sim <- simulate_redox_holobiont(
#'   n_plot = 2,
#'   n_depth = 1,
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
#' # ---- Compute cross-domain compensation index ----
#' comp <- rri_compensation_index(res)
#'
#' @export
rri_compensation_index <- function(res) {
  
  required_cols <- c("Physio", "Soil", "Micro")
  
  if (!all(required_cols %in% names(res$row_scores))) {
    stop("res$row_scores must contain Physio, Soil, and Micro columns.",
         call. = FALSE)
  }
  
  dom <- res$row_scores[, required_cols]
  dom <- stats::na.omit(dom)
  
  if (nrow(dom) < 3) {
    stop("Not enough observations to compute correlations.",
         call. = FALSE)
  }
  
  k <- ncol(dom)
  
  cor_mat <- stats::cor(dom, use = "pairwise.complete.obs")
  
  comp_val <- - (2 / (k * (k - 1))) *
    sum(cor_mat[upper.tri(cor_mat)])
  
  as.numeric(comp_val)
}