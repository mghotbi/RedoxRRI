#' @title EMS plotting theme
#'
#' @description A simple ggplot theme used for RedoxRRI visualizations.
#'
#' @param base_size Base font size
#' @return A ggplot2 theme object
#' @export
theme_ems <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold"),
      axis.title  = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid  = ggplot2::element_blank()
    )
}