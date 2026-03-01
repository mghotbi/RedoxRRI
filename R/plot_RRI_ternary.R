#' @title Plot the Holobiont Redox Resilience Ternary Diagram
#'
#' Creates a ternary plot of holobiont redox resilience based on
#' compositional domain contributions (\code{Physio}, \code{Soil}, \code{Micro})
#' and filled by per-sample \code{RRI}.
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
    tolerance = 1e-6,
    renormalize = TRUE,
    centroid_method = c("auto", "simplex_mean", "aitchison_mean")
) {
  centroid_method <- match.arg(centroid_method)
  
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
  
  # ---- ensure numeric -----------------------------------------------------
  for (nm in required) {
    ternary_df[[nm]] <- suppressWarnings(as.numeric(ternary_df[[nm]]))
  }
  
  # ---- drop undefined / degenerate rows ----------------------------------
  comps <- ternary_df[, c("Physio", "Soil", "Micro")]
  rs <- rowSums(comps)
  keep <- is.finite(rs) & rs > tolerance &
    stats::complete.cases(comps) & is.finite(ternary_df$RRI)
  ternary_df <- ternary_df[keep, , drop = FALSE]
  
  if (nrow(ternary_df) == 0) {
    stop("No valid compositional rows available for ternary plotting.")
  }
  
  # ---- closure handling ---------------------------------------------------
  comps <- ternary_df[, c("Physio", "Soil", "Micro")]
  rs2 <- rowSums(comps)
  
  if (any(abs(rs2 - 1) > tolerance)) {
    if (isTRUE(renormalize)) {
      # Renormalize (visualization-safe); guard against zeros
      rs2[rs2 == 0] <- NA_real_
      comps <- comps / rs2
      ternary_df[, c("Physio", "Soil", "Micro")] <- comps
    } else {
      stop("Physio, Soil, and Micro must sum to 1 (compositional input required).")
    }
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
    # Decide centroid method
    use_aitchison <- FALSE
    if (centroid_method == "aitchison_mean") use_aitchison <- TRUE
    if (centroid_method == "auto" && !is.null(attr(ternary_df, "clr"))) use_aitchison <- TRUE
    
    if (use_aitchison) {
      clr <- attr(ternary_df, "clr")
      
      # If attribute exists but doesn't match rows, fall back gracefully
      if (is.matrix(clr) && nrow(clr) == nrow(ternary_df) && ncol(clr) == 3) {
        # Aitchison mean: mean in clr space -> back-transform -> closure
        mu <- colMeans(clr, na.rm = TRUE)
        x <- exp(mu)
        centroid <- data.frame(
          Physio = x[1],
          Soil   = x[2],
          Micro  = x[3]
        )
        centroid <- centroid / sum(centroid)
      } else {
        # Fallback to simplex mean
        centroid <- data.frame(
          Physio = mean(ternary_df$Physio, na.rm = TRUE),
          Soil   = mean(ternary_df$Soil, na.rm = TRUE),
          Micro  = mean(ternary_df$Micro, na.rm = TRUE)
        )
        centroid <- centroid / sum(centroid)
      }
    } else {
      # Simplex arithmetic mean (legacy)
      centroid <- data.frame(
        Physio = mean(ternary_df$Physio, na.rm = TRUE),
        Soil   = mean(ternary_df$Soil, na.rm = TRUE),
        Micro  = mean(ternary_df$Micro, na.rm = TRUE)
      )
      centroid <- centroid / sum(centroid)
    }
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
  
  # ---- centroid layer -----------------------------------------------------
  if (isTRUE(show_centroid) && !is.null(centroid)) {
    p <- p +
      ggplot2::geom_point(
        data = centroid,
        ggplot2::aes(x = .data$Physio, y = .data$Soil, z = .data$Micro),
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