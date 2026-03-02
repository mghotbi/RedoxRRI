#' @title Cross-Domain Compensation Index
#'
#' @description
#' Quantifies cross-domain compensation as the negative sum of
#' off-diagonal covariances among Physio, Soil, and Micro domain scores.
#'
#' @param res An object returned by \code{rri_pipeline_st()}.
#'
#' @return A single numeric compensation index. Larger positive values
#'   indicate stronger cross-domain compensation.
#'
#' @details
#' The compensation index is defined as:
#'
#' \deqn{- \sum_{i<j} Cov(D_i, D_j)}
#'
#' where \eqn{D_i} are domain-level scores. Negative covariance
#' reflects compensatory buffering between domains.
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
  
  if (nrow(dom) < 2) {
    stop("Not enough observations to compute covariance.", call. = FALSE)
  }
  
  comp_cov <- stats::cov(dom)
  -sum(comp_cov[upper.tri(comp_cov)])
}