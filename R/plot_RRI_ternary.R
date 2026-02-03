#' Plot the Holobiont Redox Resilience Ternary Diagram
#'
#' Creates an ISME-style ternary plot of holobiont redox resilience based on
#' compositional domain contributions (\code{Physio}, \code{Soil}, \code{Micro})
#' and filled by per-sample \code{RRI}.
#'
#' @description
#' This function visualizes holobiont redox resilience in ternary space,
#' where each point represents the relative contribution of plant physiology,
#' soil redox chemistry, and microbial resilience.
#'
#' The input must contain compositional domain coordinates that sum to 1 per row.
#' Rows with undefined compositions (e.g. all-zero domains or missing values)
#' are automatically excluded.
#'
#' @param ternary_df A data frame with columns \code{Physio}, \code{Soil},
#'   \code{Micro} (compositional; sum to 1), and \code{RRI}.
#' @param point_size Numeric. Point size (default: 5).
#' @param point_alpha Numeric in (0,1]. Point transparency (default: 0.9).
#' @param palette Character. Viridis palette option
#'   (default: \code{"plasma"}).
#' @param show_subtitle Logical. Show system-level mean RRI subtitle.
#' @param show_centroid Logical. Overlay centroid (mean composition).
#' @param centroid_shape Integer. Shape for centroid marker (default: 23).
#' @param centroid_size Numeric. Size multiplier for centroid.
#' @param tolerance Numeric. Closure tolerance (default: \code{1e-6}).
#'
#' @return A \pkg{ggtern} plot object (inherits from \pkg{ggplot2}).
#'
#' @section Suggested packages:
#' \itemize{
#'   \item \pkg{ggtern}
#'   \item \pkg{ggplot2}
#'   \item \pkg{viridis}
#' }
#'
#' @examples
#' sim <- simulate_redox_holobiont(seed = 1)
#'
#' res <- rri_pipeline_st(
#'   ROS_flux = sim$ROS_flux,
#'   Eh_stability = sim$Eh_stability,
#'   micro_data = sim$micro_data,
#'   id = sim$id,
#'   direction_phys = "auto",
#'   direction_anchor_phys = "FvFm",
#'   direction_soil = "auto",
#'   direction_anchor_soil = "Eh"
#' )
#'
#' plot_RRI_ternary(res$row_scores_comp)
#'
#' @importFrom rlang .data
#' @export
plot_RRI_ternary <- function(
  ternary_df,
  point_size = 5,
  point_alpha = 0.9,
  palette = "plasma",
  show_subtitle = TRUE,
  show_centroid = TRUE,
  centroid_shape = 23,
  centroid_size = 1.4,
  tolerance = 1e-6
) {
  # ---- dependency checks --------------------------------------------------
  if (!requireNamespace("ggtern", quietly = TRUE)) {
    stop("plot_RRI_ternary() requires the {ggtern} package.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_RRI_ternary() requires the {ggplot2} package.")
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("plot_RRI_ternary() requires the {viridis} package.")
  }

  ternary_df <- as.data.frame(ternary_df)

  required <- c("Physio", "Soil", "Micro", "RRI")
  missing <- setdiff(required, names(ternary_df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  # ensure numeric
  for (nm in required) {
    ternary_df[[nm]] <- suppressWarnings(as.numeric(ternary_df[[nm]]))
  }

  # ---- drop undefined compositions ---------------------------------------
  rs <- rowSums(ternary_df[, c("Physio", "Soil", "Micro")])
  keep <- is.finite(rs) & rs > tolerance
  ternary_df <- ternary_df[keep, , drop = FALSE]

  if (nrow(ternary_df) == 0) {
    stop("No valid compositional rows available for ternary plotting.")
  }

  # ---- strict closure check ----------------------------------------------
  rs2 <- rowSums(ternary_df[, c("Physio", "Soil", "Micro")])
  if (any(abs(rs2 - 1) > tolerance)) {
    stop("Physio, Soil, and Micro must sum to 1 (compositional input required).")
  }

  # ---- subtitle -----------------------------------------------------------
  rri_index <- attr(ternary_df, "RRI_index")
  if (is.null(rri_index) || !is.finite(rri_index)) {
    rri_index <- mean(ternary_df$RRI, na.rm = TRUE)
  }

  subtitle_text <- if (isTRUE(show_subtitle)) {
    sprintf("System-level RRI index (mean): %.3f", rri_index)
  } else {
    NULL
  }

  # ---- centroid -----------------------------------------------------------
  centroid <- NULL
  if (isTRUE(show_centroid)) {
    centroid <- data.frame(
      Physio = mean(ternary_df$Physio),
      Soil   = mean(ternary_df$Soil),
      Micro  = mean(ternary_df$Micro)
    )
    centroid <- centroid / sum(centroid)
  }

  # ---- base plot ----------------------------------------------------------
  p <- ggtern::ggtern(
    data = ternary_df,
    ggplot2::aes(
      x = .data$Physio,
      y = .data$Soil,
      z = .data$Micro,
      fill = .data$RRI
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      size = point_size,
      alpha = point_alpha,
      colour = "black",
      stroke = 0.5
    ) +
    viridis::scale_fill_viridis(option = palette, direction = -1) +
    ggplot2::labs(
      title = "Holobiont Redox Resilience Triangle",
      subtitle = subtitle_text,
      x = "Plant physiology",
      y = "Soil redox chemistry",
      z = "Microbial resilience",
      fill = "RRI"
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggtern::theme_nomask() +
    ggtern::theme_showarrows() +
    ggtern::theme_clockwise()

  # ---- centroid layer (FIXED) ---------------------------------------------
  if (isTRUE(show_centroid) && !is.null(centroid)) {
    p <- p +
      ggplot2::geom_point(
        data = centroid,
        ggplot2::aes(
          x = Physio,
          y = Soil,
          z = Micro
        ),
        inherit.aes = FALSE,
        shape = centroid_shape,
        size = point_size * centroid_size,
        fill = "white",
        colour = "black",
        stroke = 1.2
      )
  }
  p
}
