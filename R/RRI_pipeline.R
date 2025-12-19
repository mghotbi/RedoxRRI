#' Holobiont Redox Resilience Index (RRI) Pipeline (Physiology + Soil + Microbial Blended Domain)
#'
#' @title
#' Redox Resilience Index (RRI) with Microbial Abundance + Network Blending
#'
#' @description
#' Computes a holobiont-level **Redox Resilience Index (RRI)** as a weighted,
#' normalized composite of three domains:
#' \enumerate{
#'   \item **Physiology** (oxidative–nitrosative buffering; \code{ROS_flux})
#'   \item **Soil** (redox stability; \code{Eh_stability})
#'   \item **Microbial resilience** (**blended** from abundance/functional structure and network organization)
#' }
#'
#' The pipeline is designed for **index construction + hypothesis testing**:
#' it converts multivariate domain measurements into **1D latent scores** (0–1),
#' applies domain-specific transforms consistent with the RRI conceptual equation,
#' then combines them with interpretable weights \eqn{w_1,w_2,w_3}.
#'
#' @details
#' ## Conceptual equation (domain-level)
#'
#' The default RRI combines normalized domain scores using:
#' \deqn{
#' \mathrm{RRI}
#' = w_1\exp(-\mathrm{ROS}^{\ast})
#' + w_2(\mathrm{Eh}^{\ast})^2
#' + w_3\mathrm{Micro}^{\ast}
#' }
#' where \eqn{(\ast)} indicates normalization to \eqn{[0,1]} and \eqn{w_1+w_2+w_3=1}.
#'
#' ### Microbial domain as a blended latent score (Option A)
#'
#' Microbial resilience is represented as a **single** domain score to preserve a
#' **3-domain ternary** (Physio/Soil/Micro) while allowing two complementary microbial aspects:
#'
#' \deqn{
#' \mathrm{Micro}^{\ast}
#' = \alpha\,\mathrm{Micro}^{\ast}_{\mathrm{abund}}
#' + (1-\alpha)\,\mathrm{Micro}^{\ast}_{\mathrm{net}}
#' }
#'
#' - \eqn{\mathrm{Micro}^{\ast}_{\mathrm{abund}}} is a latent score from \code{micro_data}
#'   (e.g., ASV/OTU abundances, MAG/pathways, functional genes).
#' - \eqn{\mathrm{Micro}^{\ast}_{\mathrm{net}}} is a latent score from \code{graph}
#'   (network organization metrics).
#' - \eqn{\alpha \in [0,1]} is controlled by \code{alpha_micro}:
#'   \itemize{
#'     \item \code{alpha_micro = 1} uses abundance/functional only
#'     \item \code{alpha_micro = 0} uses network only
#'     \item intermediate values blend both (default \code{0.5})
#'   }
#'
#' This blending is **within** the microbial domain, then the microbial domain is
#' weighted by \eqn{w_3} within the full RRI.
#'
#' ## Latent-dimension methods (recommended usage table)
#'
#' Each domain can be summarized by a **1D latent score** using different methods:
#'
#' \tabular{lll}{
#' \strong{Method} \tab \strong{Assumption / Strength} \tab \strong{Recommended when} \cr
#' \code{"mean"} \tab Equal contribution \tab Small \eqn{p}, baseline index \cr
#' \code{"scale"} \tab Standardize then average \tab Mixed units; robust baseline \cr
#' \code{"pca"} \tab Dominant linear axis \tab Moderate–strong covariance; interpretable gradient \cr
#' \code{"fa"} \tab Latent factor + noise \tab Measurement error; denoising desired \cr
#' \code{"umap"} \tab Nonlinear manifold \tab Regime shifts / thresholds / nonlinear stress responses \cr
#' \code{"nmf"} \tab Additive nonnegative structure \tab Abundances/fluxes; parts-based factors \cr
#' \code{"wgcna"} \tab Co-regulated modules \tab Soil/micro syndromes; larger \eqn{p}; correlated variables \cr
#' }
#'
#' ### Domain-specific accuracy notes (important)
#' \itemize{
#'   \item **Soil redox chemistry:** \code{"wgcna"} is often appropriate because Eh/Fe/Mn/N
#'         variables are frequently co-regulated by shared redox processes.
#'   \item **Physiology/ROS traits:** \code{"wgcna"} is appropriate only when traits reflect a
#'         coordinated stress-response syndrome; otherwise \code{"pca"} or \code{"fa"} is usually preferred.
#'   \item **Microbial abundance:** \code{"pca"}, \code{"nmf"}, \code{"wgcna"} are common choices;
#'         \code{"nmf"} is useful for nonnegative counts/abundances.
#' }
#'
#' ## Network metrics used for \eqn{\mathrm{Micro}^{\ast}_{\mathrm{net}}}
#'
#' If \code{graph} is provided, we compute:
#' \itemize{
#'   \item \eqn{Q}: Louvain modularity
#'   \item \eqn{C}: mean local clustering coefficient
#'   \item \eqn{E}: global efficiency proxy \eqn{1/\overline{d}} (inverse mean shortest path length)
#'   \item \eqn{H}: degree centralization (heterogeneity proxy)
#' }
#' These are normalized to \eqn{[0,1]} and aggregated as:
#' \deqn{
#' \mathrm{Micro}^{\ast}_{\mathrm{net}}
#' = \sqrt{\tfrac{1}{4}\left(Q^{\ast}+C^{\ast}+E^{\ast}+(1-H^{\ast})\right)}
#' }
#' (or \code{network_agg = "mean"} for a simpler arithmetic mean).
#'
#' ## Robustness and edge cases
#' \itemize{
#'   \item All per-sample outputs are scaled to \eqn{[0,1]}.
#'   \item If \code{graph} is a single \pkg{igraph} object, \eqn{\mathrm{Micro}^{\ast}_{\mathrm{net}}}
#'         is treated as **system-level** and replicated across samples.
#'   \item If \code{graph} is a list of graphs, it must have length \code{nrow(ROS_flux)}
#'         and yields sample-specific network scores.
#'   \item For \code{method = "wgcna"}: if WGCNA yields only a grey module or eigengenes cannot be computed,
#'         the method **falls back to PCA** (no hard failure).
#'   \item For \code{method = "nmf"}: data are shifted to be nonnegative if needed.
#' }
#'
#' @param ROS_flux Data frame/matrix; samples in rows, physiological/ROS traits in columns.
#' @param Eh_stability Data frame/matrix; samples in rows, soil redox chemistry variables in columns.
#' @param micro_data Optional data frame/matrix; samples in rows, microbial features in columns
#'   (e.g., ASV/OTU abundances, MAG/pathway abundances, functional genes).
#' @param graph Optional microbial network input. Either a single \pkg{igraph} object (system-level network)
#'   or a list of \pkg{igraph} objects (one per sample).
#' @param alpha_micro Numeric in \eqn{[0,1]}. Weight for microbial abundance contribution in the blended
#'   microbial score \eqn{\mathrm{Micro}^{\ast}}. Default \code{0.5}.
#' @param method_phys,method_soil,method_micro Latent extraction methods per domain. One of:
#'   \code{"pca"}, \code{"scale"}, \code{"mean"}, \code{"fa"}, \code{"umap"}, \code{"nmf"}, \code{"wgcna"}.
#' @param network_agg Character; one of \code{"equation"} (sqrt composite with \eqn{1-H}) or
#'   \code{"mean"} (mean of normalized network metrics).
#' @param w1,w2,w3 Numeric weights (>=0) that must sum to 1.
#' @param norm_method Optional legacy argument: if supplied, it overrides
#'   \code{method_phys}, \code{method_soil}, and \code{method_micro}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{Physio}{Compositional physiological component (0–1).}
#'   \item{Soil}{Compositional soil component (0–1).}
#'   \item{Micro}{Compositional blended microbial component (0–1).}
#'   \item{RRI}{Per-sample Redox Resilience Index (0–1).}
#'   \item{RRI_index}{System-level mean RRI across samples.}
#'   \item{Micro_abundance}{Unblended microbial abundance latent score (0–1), if \code{micro_data} provided.}
#'   \item{Micro_network}{Unblended microbial network score (0–1), if \code{graph} provided.}
#' }
#' The mean RRI is also stored as \code{attr(df, "RRI_index")}.
#'
#' @section Package dependency guidance:
#' \itemize{
#'   \item **Imports:** \pkg{stats}, \pkg{igraph}
#'   \item **Suggests:** \pkg{WGCNA}, \pkg{psych}, \pkg{uwot}, \pkg{NMF}
#' }
#' Optional methods are activated only when those packages are installed.
#'
#' @examples
#' # ---- Physiology + soil + microbial abundance + network ------------------
#' set.seed(1)
#' ROS_flux <- data.frame(
#'   SPAD    = rnorm(30, 35, 3),
#'   FvFm    = rnorm(30, 0.78, 0.03),
#'   PhiPSII = rnorm(30, 0.40, 0.05),
#'   NPQ     = rnorm(30, 1.2, 0.2)
#' )
#'
#' Eh_stability <- data.frame(
#'   Eh       = runif(30, 0.3, 0.9),
#'   Fe2.Fe3  = rbeta(30, 3, 5),
#'   Mn2.Mn4  = rbeta(30, 2, 6),
#'   NH4.NO3  = rnorm(30, 1.5, 0.4)
#' )
#'
#' # microbial abundance/features (e.g., ASVs)
#' micro_data <- matrix(rexp(30 * 50, rate = 2), nrow = 30)
#' colnames(micro_data) <- paste0("ASV", seq_len(ncol(micro_data)))
#'
#' # microbial network (system-level)
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   g <- igraph::sample_smallworld(1, 40, 5, 0.10)
#'   g <- igraph::simplify(g)
#'
#'   out <- RRI_pipeline(
#'     ROS_flux,
#'     Eh_stability,
#'     micro_data   = micro_data,
#'     graph        = g,
#'     alpha_micro  = 0.6,         # 60% abundance, 40% network
#'     method_phys  = "pca",
#'     method_soil  = "wgcna",
#'     method_micro = "pca",
#'     network_agg  = "equation",
#'     w1 = 0.4, w2 = 0.35, w3 = 0.25
#'   )
#'   head(out)
#'   attr(out, "RRI_index")
#' }
#'
#' @importFrom stats prcomp
#' @importFrom igraph simplify cluster_louvain modularity transitivity distances centr_degree
#' @export
RRI_pipeline <- function(
    ROS_flux,
    Eh_stability,
    micro_data = NULL,
    graph = NULL,
    alpha_micro = 0.5,
    method_phys = "pca",
    method_soil = "pca",
    method_micro = "pca",
    network_agg = c("equation", "mean"),
    w1 = 0.4,
    w2 = 0.35,
    w3 = 0.25,
    norm_method = NULL
) {
  network_agg <- match.arg(network_agg)

  # legacy override: one method for all latent domains
  if (!is.null(norm_method)) {
    method_phys  <- norm_method
    method_soil  <- norm_method
    method_micro <- norm_method
  }

  # ---- helpers -------------------------------------------------------------
  .as_numeric_df <- function(df) {
    df <- as.data.frame(df)
    df[] <- lapply(df, function(x) suppressWarnings(as.numeric(as.character(x))))
    df
  }

  .scale_01 <- function(v) {
    if (all(is.na(v))) return(rep(NA_real_, length(v)))
    rng <- range(v, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0.5, length(v)))
    (v - rng[1]) / diff(rng)
  }

  .check_weights <- function(w1, w2, w3) {
    w <- c(w1, w2, w3)
    if (any(!is.finite(w)) || any(w < 0)) stop("Weights must be finite and >= 0.")
    if (abs(sum(w) - 1) > 1e-8) stop("Weights must sum to 1. Got: ", sum(w))
    invisible(TRUE)
  }

  .align_rows <- function(a, b, c = NULL) {
    a <- .as_numeric_df(a)
    b <- .as_numeric_df(b)

    if (nrow(a) != nrow(b)) {
      stop("ROS_flux and Eh_stability must have the same number of rows (samples).")
    }

    if (!is.null(c)) {
      c <- .as_numeric_df(c)
      if (nrow(c) != nrow(a)) {
        stop("micro_data must have the same number of rows as ROS_flux.")
      }
    }

    list(a = a, b = b, c = c)
  }

  latent_dimension <- function(df, method = "pca") {
    df <- .as_numeric_df(df)

    if (ncol(df) == 1L) return(.scale_01(df[[1L]]))

    if (method == "scale") {
      v <- rowMeans(scale(df), na.rm = TRUE)

    } else if (method == "mean") {
      v <- rowMeans(df, na.rm = TRUE)

    } else if (method == "pca") {
      v <- stats::prcomp(df, center = TRUE, scale. = TRUE)$x[, 1]

    } else if (method == "fa") {
      if (!requireNamespace("psych", quietly = TRUE)) {
        stop("method = 'fa' requires the 'psych' package.")
      }
      fa <- psych::fa(df, nfactors = 1, rotate = "none", scores = "regression")
      v <- as.numeric(fa$scores[, 1])

    } else if (method == "umap") {
      if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("method = 'umap' requires the 'uwot' package.")
      }
      v <- as.numeric(uwot::umap(df, n_components = 1))

    } else if (method == "nmf") {
      if (!requireNamespace("NMF", quietly = TRUE)) {
        stop("method = 'nmf' requires the 'NMF' package.")
      }
      x <- as.matrix(df)

      # shift to nonnegative
      mins  <- apply(x, 2, min, na.rm = TRUE)
      shift <- pmax(0, -mins)
      x <- sweep(x, 2, shift, "+")

      nm <- NMF::nmf(x, rank = 1, .options = "N")
      v <- as.numeric(t(NMF::coef(nm)))

    } else if (method == "wgcna") {
      if (!requireNamespace("WGCNA", quietly = TRUE)) {
        stop("method = 'wgcna' requires the 'WGCNA' package.")
      }

      df2 <- df
      gsg <- WGCNA::goodSamplesGenes(df2, verbose = 0)
      if (!gsg$allOK) df2 <- df2[gsg$goodSamples, gsg$goodGenes, drop = FALSE]

      # small/noisy -> PCA fallback
      if (nrow(df2) < 4L || ncol(df2) < 3L) {
        return(.scale_01(stats::prcomp(df, center = TRUE, scale. = TRUE)$x[, 1]))
      }

      sft <- suppressWarnings(WGCNA::pickSoftThreshold(df2, verbose = 0))
      power <- sft$powerEstimate
      if (is.na(power) || !is.finite(power)) power <- 6

      net <- suppressWarnings(
        WGCNA::blockwiseModules(
          df2,
          power          = power,
          TOMType        = "unsigned",
          minModuleSize  = 5,
          mergeCutHeight = 0.25,
          numericLabels  = TRUE,
          verbose        = 0
        )
      )

      cols <- net$colors
      if (length(unique(cols)) <= 1L) {
        # only grey -> PCA fallback on df2
        v <- stats::prcomp(df2, center = TRUE, scale. = TRUE)$x[, 1]
        out <- rep(NA_real_, nrow(df))
        out[as.integer(rownames(df2))] <- v
        return(.scale_01(out))
      }

      eig <- suppressWarnings(WGCNA::moduleEigengenes(df2, cols)$eigengenes)
      if (is.null(eig) || ncol(eig) < 1L) {
        v <- stats::prcomp(df2, center = TRUE, scale. = TRUE)$x[, 1]
        out <- rep(NA_real_, nrow(df))
        out[as.integer(rownames(df2))] <- v
        return(.scale_01(out))
      }

      v <- as.numeric(eig[, 1])
      out <- rep(NA_real_, nrow(df))
      out[as.integer(rownames(df2))] <- v
      return(.scale_01(out))

    } else {
      stop("Unknown latent method: ", method)
    }

    .scale_01(v)
  }

  network_scalar <- function(g, agg = "equation") {
    g <- igraph::simplify(g)

    comm <- igraph::cluster_louvain(g)
    Q <- igraph::modularity(comm)

    C_vec <- igraph::transitivity(g, type = "local", isolates = "zero")
    C <- mean(C_vec, na.rm = TRUE)

    D <- igraph::distances(g)
    D[!is.finite(D)] <- NA_real_
    Eglob <- 1 / mean(D, na.rm = TRUE)

    # degree centralization as heterogeneity proxy
    H <- igraph::centralization.degree(g)$centralization

    # normalize (safe bounds for Q,C,H; scale for Eglob)
    Qs <- pmin(pmax(Q, 0), 1)
    Cs <- pmin(pmax(C, 0), 1)
    Hs <- pmin(pmax(H, 0), 1)
    Es <- .scale_01(Eglob)

    metrics <- c(Qs, Cs, Es, (1 - Hs))
    if (agg == "mean") return(mean(metrics, na.rm = TRUE))
    sqrt(mean(metrics, na.rm = TRUE))
  }

  # ---- start ---------------------------------------------------------------
  .check_weights(w1, w2, w3)

  if (!is.finite(alpha_micro) || alpha_micro < 0 || alpha_micro > 1) {
    stop("alpha_micro must be a finite value between 0 and 1.")
  }

  aligned <- .align_rows(ROS_flux, Eh_stability, micro_data)
  ROS_flux <- aligned$a
  Eh_stability <- aligned$b
  micro_data <- aligned$c
  n <- nrow(ROS_flux)

  # --- domain latents (0..1) -----------------------------------------------
  ros_lat <- latent_dimension(ROS_flux, method = method_phys)
  eh_lat  <- latent_dimension(Eh_stability, method = method_soil)

  # --- microbial abundance latent (optional) -------------------------------
  micro_abund <- rep(NA_real_, n)
  if (!is.null(micro_data)) {
    micro_abund <- latent_dimension(micro_data, method = method_micro)
  }

  # --- microbial network latent (optional) ---------------------------------
  micro_net <- rep(NA_real_, n)
  if (!is.null(graph)) {
    if (inherits(graph, "igraph")) {
      micro_net <- rep(network_scalar(graph, agg = network_agg), n)
    } else if (is.list(graph)) {
      if (length(graph) != n) stop("If graph is a list, it must have length nrow(ROS_flux).")
      micro_net <- vapply(graph, function(g) network_scalar(g, agg = network_agg), numeric(1))
    } else {
      stop("graph must be NULL, an igraph object, or a list of igraph objects.")
    }
    micro_net <- .scale_01(micro_net)
  }

  # require at least one microbial source
  if (all(is.na(micro_abund)) && all(is.na(micro_net))) {
    stop("Provide microbial input via micro_data and/or graph to compute the microbial domain.")
  }

  # --- blended microbial domain (Micro*) -----------------------------------
  # If only one microbial source is provided, Micro* becomes that source.
  if (all(is.na(micro_abund))) {
    micro_lat <- micro_net
  } else if (all(is.na(micro_net))) {
    micro_lat <- micro_abund
  } else {
    micro_lat <- .scale_01(alpha_micro * micro_abund + (1 - alpha_micro) * micro_net)
  }

  # --- apply RRI transforms and combine ------------------------------------
  term_ros <- exp(-ros_lat)
  term_eh  <- eh_lat^2

  rri <- .scale_01(w1 * term_ros + w2 * term_eh + w3 * micro_lat)

  # --- assemble output ------------------------------------------------------
  out <- data.frame(
    Physio = ros_lat,
    Soil   = eh_lat,
    Micro  = micro_lat,
    RRI    = rri,
    Micro_abundance = micro_abund,
    Micro_network   = micro_net
  )

  # compositional projection for ternary visualization uses blended Micro
  comp <- out[, c("Physio", "Soil", "Micro")]
  comp_sum <- rowSums(comp)
  comp_sum[comp_sum == 0] <- NA_real_
  out[, c("Physio", "Soil", "Micro")] <- comp / comp_sum

  rri_index <- mean(out$RRI, na.rm = TRUE)
  out$RRI_index <- rri_index
  attr(out, "RRI_index") <- rri_index
  attr(out, "alpha_micro") <- alpha_micro

  out
}

# library(igraph)
# set.seed(1)
#
# n <- 30
#
# # --- Physiology / ROS domain (samples x variables) -------------------------
# ROS_flux <- data.frame(
#   SPAD    = rnorm(n, 35, 3),
#   FvFm    = rnorm(n, 0.78, 0.03),
#   PhiPSII = rnorm(n, 0.40, 0.05),
#   NPQ     = rnorm(n, 1.2, 0.2)
# )
#
# # --- Soil redox chemistry domain ------------------------------------------
# Eh_stability <- data.frame(
#   Eh       = runif(n, 0.3, 0.9),
#   Fe2.Fe3  = rbeta(n, 3, 5),
#   Mn2.Mn4  = rbeta(n, 2, 6),
#   NH4.NO3  = rnorm(n, 1.5, 0.4)
# )
#
# # --- Microbial abundance/features (nonnegative) ----------------------------
# # Example: 60 ASVs/genes; zero-inflated-ish counts
# p_micro <- 60
# micro_data <- matrix(rpois(n * p_micro, lambda = 5), nrow = n, ncol = p_micro)
# micro_data[sample(length(micro_data), size = round(0.25 * length(micro_data)))] <- 0
# micro_data <- as.data.frame(micro_data)
# colnames(micro_data) <- paste0("ASV", seq_len(p_micro))
#
# # --- Microbial network (system-level graph) --------------------------------
# g <- igraph::sample_smallworld(1, 40, 5, 0.10)
# g <- igraph::simplify(g)
# tdf <- RRI_pipeline(
#   ROS_flux,
#   Eh_stability,
#   micro_data   = micro_data,  # abundance/features
#   graph        = g,           # network
#   alpha_micro  = 0.6,         # 60% abundance, 40% network
#   method_micro = "pca",
#   method_soil  = "wgcna"
# )
#
# attr(tdf, "RRI_index")
# head(tdf)
# summary(tdf$RRI)
# set.seed(1)
# g_list <- replicate(
#   nrow(ROS_flux),
#   {
#     gg <- igraph::rewire(g, with = igraph::keeping_degseq(niter = 50))
#     igraph::simplify(gg)
#   },
#   simplify = FALSE
# )
#
# tdf2 <- RRI_pipeline(
#   ROS_flux,
#   Eh_stability,
#   micro_data  = micro_data,
#   graph       = g_list,
#   alpha_micro = 0.6,
#   method_soil = "wgcna"
# )
#
# sd(tdf2$Micro_network)
# cor(tdf2$Micro_abundance, tdf2$Micro_network, use = "pairwise.complete.obs")

