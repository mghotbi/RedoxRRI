#' @title Latent Recovery Correlation for RedoxRRI
#'
#' @description
#' Computes the Pearson correlation between the per-sample Redox Resilience
#' Index (RRI) and a known latent truth vector. This function is intended
#' for simulation-based validation only and should not be interpreted as
#' predictive accuracy.
#'
#' @param res An object returned by \code{rri_pipeline_st()}.
#' @param latent_truth Numeric vector of true latent values.
#' @param method Correlation method. One of \code{"pearson"},
#'   \code{"spearman"}, or \code{"kendall"}.
#'
#' @return A single numeric correlation coefficient.
#'
#' @details
#' This function is designed for simulation benchmarking. In empirical
#' datasets, no latent truth exists and this metric should not be used.
#'
#' @export
rri_latent_correlation <- function(res,
                                   latent_truth,
                                   method = c("pearson", "spearman", "kendall")) {
  
  method <- match.arg(method)
  
  if (is.null(res$row_scores$RRI)) {
    stop("res must contain a row_scores$RRI column.", call. = FALSE)
  }
  
  rri <- as.numeric(res$row_scores$RRI)
  latent_truth <- as.numeric(latent_truth)
  
  if (length(rri) != length(latent_truth)) {
    stop("RRI and latent_truth must have equal length.", call. = FALSE)
  }
  
  stats::cor(rri, latent_truth, method = method, use = "complete.obs")
}