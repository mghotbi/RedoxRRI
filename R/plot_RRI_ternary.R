#' Plot the Holobiont Redox Resilience Ternary Diagram
#'
#' @title
#' Holobiont Redox Resilience Triangle
#'
#' @description
#' Creates an ISME-style ternary plot of holobiont redox resilience based on the
#' output of \code{\link{RRI_pipeline}}. Points are positioned by the compositional
#' domain contributions (\code{Physio}, \code{Soil}, \code{Micro}) and colored by the
#' per-sample \code{RRI}. The system-level mean RRI (\code{attr(x, "RRI_index")}) is
#' shown in the subtitle.
#'
#' @details
#' This function expects a data frame with columns:
#' \itemize{
#'   \item \code{Physio}, \code{Soil}, \code{Micro} (compositional; should sum to 1 per row),
#'   \item \code{RRI} (0–1).
#' }
#'
#' Plotting packages are treated as **optional** dependencies (listed in \code{Suggests}).
#' If not installed, the function stops with a clear message.
#'
#' @param ternary_df A data frame produced by \code{\link{RRI_pipeline}}.
#' @param point_size Numeric; point size (default \code{5}).
#' @param point_alpha Numeric; point alpha (default \code{0.90}).
#' @param palette Character; viridis palette option (default \code{"plasma"}).
#' @param show_subtitle Logical; include system mean RRI subtitle (default \code{TRUE}).
#'
#' @return A \code{ggtern} object.
#'
#' @section Suggested packages:
#' \itemize{
#'   \item \pkg{ggtern}
#'   \item \pkg{ggplot2}
#'   \item \pkg{viridis}
#' }
#'
#' @examples
#' # Requires ggtern + ggplot2 + viridis installed
#' # out <- RRI_pipeline(...)
#' # p <- plot_RRI_ternary(out)
#' # print(p)
#' @importFrom rlang .data
#' @export
plot_RRI_ternary <- function(
    ternary_df,
    point_size = 5,
    point_alpha = 0.90,
    palette = "plasma",
    show_subtitle = TRUE
) {
  if (!requireNamespace("ggtern", quietly = TRUE)) {
    stop("plot_RRI_ternary() requires the 'ggtern' package (Suggested). Install it to plot.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_RRI_ternary() requires the 'ggplot2' package (Suggested). Install it to plot.")
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("plot_RRI_ternary() requires the 'viridis' package (Suggested). Install it to plot.")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("plot_RRI_ternary() requires the 'grid' package (part of base R distributions).")
  }

  required_cols <- c("Physio", "Soil", "Micro", "RRI")
  missing_cols <- setdiff(required_cols, colnames(ternary_df))
  if (length(missing_cols) > 0) {
    stop(
      "ternary_df is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      ". Input should be the output of RRI_pipeline()."
    )
  }

  rri_index <- attr(ternary_df, "RRI_index")
  if (is.null(rri_index)) {
    rri_index <- mean(ternary_df$RRI, na.rm = TRUE)
  }

  subtitle_text <- if (isTRUE(show_subtitle) && is.finite(rri_index)) {
    sprintf("System-level RRI index (mean): %.3f", rri_index)
  } else {
    NULL
  }

  ggtern::ggtern(
    data = ternary_df,
    ggplot2::aes(x = Physio, y = Soil, z = Micro, colour = RRI)
  ) +
    # outer stroke ring for crisp points
    ggplot2::geom_point(
      size   = point_size,
      alpha  = point_alpha * 0.75,
      shape  = 21,
      stroke = 0.6,
      colour = "black",
      fill   = NA
    ) +
    # colored fill layer
    ggplot2::geom_point(
      size  = point_size,
      alpha = point_alpha
    ) +
    viridis::scale_color_viridis(
      option    = palette,
      direction = -1
    ) +
    ggplot2::labs(
      title    = "Holobiont Redox Resilience Triangle",
      subtitle = subtitle_text,
      x = "Plant physiology",
      y = "Soil redox chemistry",
      z = "Microbial resilience",
      colour = "RRI"
    ) +
    ggplot2::theme_bw(base_size = 15) +
    ggtern::theme_nomask() +
    ggtern::theme_showarrows() +
    ggtern::theme_clockwise() +
    ggplot2::theme(
      # Remove ggplot2's stray cartesian axis label that can appear as "y"
      axis.title.y = ggplot2::element_blank(),

      # Optional: remove redundant top tern title (ggtern can duplicate)
      tern.axis.title.T = ggplot2::element_blank(),

      # keep bottom + right at reasonable positions
      tern.axis.title.L = ggplot2::element_text(hjust = 1.1),
      tern.axis.title.R = ggplot2::element_text(hjust = -0.1),

      plot.margin = grid::unit(c(1.4, 1.2, 1.2, 1.2), "lines"),

      plot.title = ggplot2::element_text(
        face  = "bold",
        size  = 18,
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size  = 12,
        hjust = 0.5
      ),

      axis.title = ggplot2::element_text(face = "bold"),
      axis.text  = ggplot2::element_text(size = 11),

      legend.position = "right",
      legend.title    = ggplot2::element_text(face = "bold", size = 12),
      legend.text     = ggplot2::element_text(size = 10)
    )
}



