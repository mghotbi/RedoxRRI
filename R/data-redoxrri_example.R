#' Example data for RedoxRRI
#'
#' A small simulated dataset bundled with \pkg{RedoxRRI} to demonstrate the
#' full workflow: computing the Redox Resilience Index (RRI) and creating a
#' ternary visualization.
#'
#' @format A named list with components:
#' \describe{
#'   \item{ROS_flux}{data.frame; samples (rows) by physiological/ROS traits (columns).}
#'   \item{Eh_stability}{data.frame; samples (rows) by soil redox chemistry variables (columns).}
#'   \item{micro_data}{matrix/data.frame; samples (rows) by microbial features (columns).}
#'   \item{graph}{An \pkg{igraph} object representing a microbial network.}
#' }
#'
#' @source Simulated data generated for examples and unit tests.
#'
#' @examples
#' data("redoxrri_example")
#' str(redoxrri_example)
#'
#' out <- RRI_pipeline(
#'   ROS_flux     = redoxrri_example$ROS_flux,
#'   Eh_stability = redoxrri_example$Eh_stability,
#'   micro_data   = redoxrri_example$micro_data,
#'   graph        = redoxrri_example$graph,
#'   alpha_micro  = 0.6
#' )
#'
#' head(out)
#' attr(out, "RRI_index")
#'
#' @docType data
#' @name redoxrri_example
NULL

"redoxrri_example"
