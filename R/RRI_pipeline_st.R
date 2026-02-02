#' Holobiont Redox Resilience Index (RRI) with spatio-temporal dynamics
#'
#' Computes a holobiont-level Redox Resilience Index (RRI) by integrating
#' plant physiology, soil redox chemistry, and microbial resilience into a unified,
#' directionally identifiable index with explicit support for spatial and temporal designs.
#'
#' @description
#' This function constructs a holobiont-scale resilience index from three domains:
#' \enumerate{
#'   \item \strong{Physiology} (plant oxidative buffering and stress-response traits)
#'   \item \strong{Soil} (redox chemistry and stability proxies)
#'   \item \strong{Microbial} resilience, optionally \strong{blended} from:
#'     \itemize{
#'       \item microbial abundance or functional features (\code{micro_data})
#'       \item microbial network organization (\code{graph})
#'     }
#' }
#'
#' @details
#' ## Returned domain representations
#' The output includes both:
#' \itemize{
#'   \item \strong{Absolute domain scores} (\code{Physio}, \code{Soil}, \code{Micro}) scaled to \eqn{[0,1]}.
#'   \item \strong{Compositional domain contributions} (\code{Physio}, \code{Soil}, \code{Micro})
#'   projected to sum to 1 (ternary-ready).
#' }
#'
#' ## Direction identifiability
#' Latent methods (e.g., PCA/FA) have arbitrary sign. When \code{direction_* = "auto"},
#' the latent score is flipped (if needed) so that higher latent values align with higher
#' values of the specified anchor variable.
#'
#' ## Missing data handling
#' Robust median imputation is applied \emph{within strata} (defined by \code{scale_by})
#' before latent extraction. This avoids failures in SVD-based methods and reduces bias
#' under MNAR sensor dropout and redox shocks.
#'
#' ## Spatio-temporal modes
#' \describe{
#'   \item{\code{"snapshot"}}{Independent per-sample RRI (default).}
#'   \item{\code{"rolling"}}{Rolling-window dynamic RRI within \code{group_cols} ordered by \code{time_col}.}
#'   \item{\code{"event"}}{Baseline–pulse–recovery summaries by group using \code{event_col}.}
#' }
#'
#' @param ROS_flux Data frame or matrix; samples in rows, physiological traits in columns.
#' @param Eh_stability Data frame or matrix; samples in rows, soil redox variables in columns.
#' @param micro_data Optional data frame or matrix; samples in rows, microbial features in columns.
#' @param graph Optional microbial network input: either a single \pkg{igraph} object (system-level)
#'   or a list of \pkg{igraph} objects (one per sample).
#'
#' @param id Optional data frame describing sample identities. If provided, must have
#'   \code{nrow(id) == nrow(ROS_flux)}. Typical columns: \code{plot}, \code{depth}, \code{plant_id}, \code{time}.
#' @param time_col Character. Name of time column in \code{id}. Required for \code{mode != "snapshot"}.
#' @param group_cols Character vector. Columns in \code{id} defining spatial/experimental units.
#'
#' @param mode One of \code{"snapshot"}, \code{"rolling"}, \code{"event"}.
#' @param window Integer >= 2. Rolling window size for \code{mode = "rolling"}.
#' @param align Character. Rolling alignment: \code{"right"}, \code{"center"}, \code{"left"}.
#'
#' @param event_col Character. Column in \code{id} containing event labels (for \code{mode = "event"}).
#' @param baseline_label Character. Label defining baseline period.
#' @param recovery_labels Character vector defining recovery periods.
#'
#' @param alpha_micro Numeric in \eqn{[0,1]}. Blend weight:
#'   \deqn{\mathrm{Micro} = \alpha \cdot \mathrm{abundance} + (1-\alpha) \cdot \mathrm{network}.}
#'
#' @param method_phys,method_soil,method_micro Latent extraction methods:
#'   \code{"pca"}, \code{"scale"}, \code{"mean"}, \code{"fa"}, \code{"umap"}, \code{"nmf"}, \code{"wgcna"}.
#'
#' @param direction_phys,direction_soil,direction_micro Direction semantics:
#'   \code{"higher_is_better"}, \code{"lower_is_better"}, \code{"auto"}.
#' @param direction_anchor_phys,direction_anchor_soil,direction_anchor_micro
#'   Optional anchor column names used when corresponding direction is \code{"auto"}.
#'
#' @param scale_by Character vector of columns in \code{id} within which to scale domain scores to \eqn{[0,1]}.
#'   Defaults to \code{group_cols}. Set \code{NULL} for global scaling.
#'
#' @param network_agg Aggregation for network metrics: \code{"equation"} or \code{"mean"}.
#'
#' @param w1,w2,w3 Nonnegative domain weights that must sum to 1.
#'
#' @param add_coupling Logical. If TRUE, adds a coupling/coherence term among domains.
#' @param coupling_weight Numeric in \eqn{[0,1]}. Weight for coupling term; domain weights are renormalized
#'   to sum to \eqn{1 - \mathrm{coupling\_weight}}.
#' @param coupling_fun Coupling function: \code{"geometric_mean"} or \code{"agreement"}.
#'
#' @param norm_method Optional legacy override; if supplied, overrides \code{method_*}.
#'
#' @return An object of class \code{"RRI"} (a named list) with:
#' \describe{
#'   \item{row_scores}{Absolute domain scores and \code{RRI} (all scaled to \eqn{[0,1]}).}
#'   \item{row_scores_comp}{Compositional domain contributions and \code{RRI} (Physio + Soil + Micro = 1 when defined).}
#'   \item{dyn_scores}{Dynamic summaries for rolling/event modes; \code{NULL} for snapshot mode.}
#'   \item{meta}{Settings and system-level mean RRI (\code{rri_index}).}
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
#'   group_cols = c("plot", "depth"),
#'   scale_by = c("plot", "depth"),
#'   direction_phys = "auto",
#'   direction_anchor_phys = "FvFm",
#'   direction_soil = "auto",
#'   direction_anchor_soil = "Eh",
#'   direction_micro = "higher_is_better"
#' )
#'
#' head(res$row_scores)
#' head(res$row_scores_comp)
#'
#' @importFrom stats prcomp lm sd coef cor median
#' @export
rri_pipeline_st <- function(
  ROS_flux,
  Eh_stability,
  micro_data = NULL,
  graph = NULL,
  id = NULL,
  time_col = NULL,
  group_cols = NULL,
  mode = c("snapshot", "rolling", "event"),
  window = 3,
  align = c("right", "center", "left"),
  event_col = NULL,
  baseline_label = "pre",
  recovery_labels = "recovery",
  alpha_micro = 0.5,
  method_phys = "pca",
  method_soil = "pca",
  method_micro = "pca",
  direction_phys = c("auto", "higher_is_better", "lower_is_better"),
  direction_soil = c("auto", "higher_is_better", "lower_is_better"),
  direction_micro = c("auto", "higher_is_better", "lower_is_better"),
  direction_anchor_phys = NULL,
  direction_anchor_soil = NULL,
  direction_anchor_micro = NULL,
  scale_by = NULL,
  network_agg = c("equation", "mean"),
  w1 = 0.4,
  w2 = 0.35,
  w3 = 0.25,
  add_coupling = FALSE,
  coupling_weight = 0,
  coupling_fun = c("geometric_mean", "agreement"),
  norm_method = NULL
) {
  mode <- match.arg(mode)
  align <- match.arg(align)
  network_agg <- match.arg(network_agg)
  coupling_fun <- match.arg(coupling_fun)
  direction_phys <- match.arg(direction_phys)
  direction_soil <- match.arg(direction_soil)
  direction_micro <- match.arg(direction_micro)

  if (!is.null(norm_method)) {
    method_phys <- norm_method
    method_soil <- norm_method
    method_micro <- norm_method
  }

  # ---- validations --------------------------------------------------------
  if (missing(ROS_flux) || missing(Eh_stability)) {
    stop("Provide `ROS_flux` and `Eh_stability`.")
  }

  ROS_flux <- as.data.frame(ROS_flux)
  Eh_stability <- as.data.frame(Eh_stability)

  if (nrow(ROS_flux) != nrow(Eh_stability)) {
    stop("ROS_flux and Eh_stability must have the same number of rows (samples).")
  }
  n <- nrow(ROS_flux)

  if (!is.null(id)) {
    id <- as.data.frame(id)
    if (nrow(id) != n) stop("`id` must have the same number of rows as `ROS_flux`.")
  }
  if (!is.null(micro_data)) {
    micro_data <- as.data.frame(micro_data)
    if (nrow(micro_data) != n) stop("`micro_data` must have the same number of rows as `ROS_flux`.")
  }

  if (any(!is.finite(c(w1, w2, w3))) || any(c(w1, w2, w3) < 0)) {
    stop("w1, w2, w3 must be finite and >= 0.")
  }
  if (abs(w1 + w2 + w3 - 1) > 1e-8) stop("w1 + w2 + w3 must sum to 1.")

  if (!is.finite(alpha_micro) || alpha_micro < 0 || alpha_micro > 1) stop("alpha_micro must be in [0,1].")
  if (!is.finite(coupling_weight) || coupling_weight < 0 || coupling_weight > 1) stop("coupling_weight must be in [0,1].")

  if (mode != "snapshot") {
    if (is.null(id) || is.null(time_col) || !time_col %in% names(id)) {
      stop("For mode != 'snapshot', provide `id` with a valid `time_col`.")
    }
    if (is.null(group_cols) || any(!group_cols %in% names(id))) {
      stop("For mode != 'snapshot', provide valid `group_cols` in `id`.")
    }
    if (!is.numeric(window) || length(window) != 1 || window < 2) stop("`window` must be a single integer >= 2.")
  }

  if (is.null(scale_by) && !is.null(group_cols)) scale_by <- group_cols
  if (!is.null(scale_by) && (is.null(id) || any(!scale_by %in% names(id)))) {
    stop("If `scale_by` is provided, `id` must be provided and contain all `scale_by` columns.")
  }

  # ---- helpers ------------------------------------------------------------
  as_numeric_df <- function(x) {
    x <- as.data.frame(x)
    x[] <- lapply(x, function(v) suppressWarnings(as.numeric(as.character(v))))
    x
  }

  impute_median_by <- function(df, id, by) {
    df <- as_numeric_df(df)

    if (is.null(by) || is.null(id) || length(by) == 0) {
      for (j in seq_len(ncol(df))) {
        v <- df[[j]]
        if (all(is.na(v))) next
        v[is.na(v)] <- stats::median(v, na.rm = TRUE)
        df[[j]] <- v
      }
      return(df)
    }

    key <- interaction(id[, by, drop = FALSE], drop = TRUE, lex.order = TRUE)
    for (k in unique(key)) {
      idx <- which(key == k)
      for (j in seq_len(ncol(df))) {
        v <- df[[j]][idx]
        if (all(is.na(v))) next
        v[is.na(v)] <- stats::median(v, na.rm = TRUE)
        df[[j]][idx] <- v
      }
    }

    for (j in seq_len(ncol(df))) {
      v <- df[[j]]
      if (anyNA(v) && !all(is.na(v))) {
        v[is.na(v)] <- stats::median(v, na.rm = TRUE)
        df[[j]] <- v
      }
    }

    df
  }

  scale_01 <- function(v) {
    if (all(is.na(v))) {
      return(rep(NA_real_, length(v)))
    }
    r <- range(v, na.rm = TRUE)
    if (!all(is.finite(r)) || diff(r) == 0) {
      return(rep(0.5, length(v)))
    }
    (v - r[1]) / diff(r)
  }

  scale_01_by <- function(v, id, by) {
    if (is.null(by) || is.null(id) || length(by) == 0) {
      return(scale_01(v))
    }
    key <- interaction(id[, by, drop = FALSE], drop = TRUE, lex.order = TRUE)
    out <- v
    for (k in unique(key)) {
      idx <- which(key == k)
      out[idx] <- scale_01(v[idx])
    }
    out
  }

  orient_latent <- function(latent, df, direction, anchor) {
    if (direction == "higher_is_better") {
      return(latent)
    }
    if (direction == "lower_is_better") {
      return(-latent)
    }
    if (is.null(anchor) || !anchor %in% names(df)) {
      return(latent)
    }
    a <- df[[anchor]]
    ok <- is.finite(latent) & is.finite(a)
    if (!any(ok)) {
      return(latent)
    }
    if (stats::cor(latent[ok], a[ok]) < 0) -latent else latent
  }

  latent_dimension <- function(df, method) {
    df <- as_numeric_df(df)
    if (ncol(df) == 1L) {
      return(df[[1L]])
    }

    if (method == "mean") {
      return(rowMeans(df, na.rm = TRUE))
    }
    if (method == "scale") {
      return(rowMeans(scale(df), na.rm = TRUE))
    }
    if (method == "pca") {
      return(stats::prcomp(df, center = TRUE, scale. = TRUE)$x[, 1])
    }

    if (method == "fa") {
      if (!requireNamespace("psych", quietly = TRUE)) stop("method = 'fa' requires {psych}.")
      fa <- psych::fa(df, nfactors = 1, rotate = "none", scores = "regression")
      return(as.numeric(fa$scores[, 1]))
    }

    if (method == "umap") {
      if (!requireNamespace("uwot", quietly = TRUE)) stop("method = 'umap' requires {uwot}.")
      return(as.numeric(uwot::umap(df, n_components = 1)))
    }

    if (method == "nmf") {
      if (!requireNamespace("NMF", quietly = TRUE)) stop("method = 'nmf' requires {NMF}.")
      x <- as.matrix(df)
      mins <- apply(x, 2, min, na.rm = TRUE)
      shift <- pmax(0, -mins)
      x <- sweep(x, 2, shift, "+")
      nm <- NMF::nmf(x, rank = 1, .options = "N")
      return(as.numeric(t(NMF::coef(nm))))
    }

    if (method == "wgcna") {
      if (!requireNamespace("WGCNA", quietly = TRUE)) stop("method = 'wgcna' requires {WGCNA}.")
      gsg <- WGCNA::goodSamplesGenes(df, verbose = 0)
      df2 <- df
      if (!gsg$allOK) df2 <- df2[gsg$goodSamples, gsg$goodGenes, drop = FALSE]

      if (nrow(df2) < 4L || ncol(df2) < 3L) {
        return(stats::prcomp(df, center = TRUE, scale. = TRUE)$x[, 1])
      }

      sft <- suppressWarnings(WGCNA::pickSoftThreshold(df2, verbose = 0))
      power <- sft$powerEstimate
      if (is.na(power) || !is.finite(power)) power <- 6

      net <- suppressWarnings(
        WGCNA::blockwiseModules(
          df2,
          power = power,
          TOMType = "unsigned",
          minModuleSize = 5,
          mergeCutHeight = 0.25,
          numericLabels = TRUE,
          verbose = 0
        )
      )

      cols <- net$colors
      if (length(unique(cols)) <= 1L) {
        out <- rep(NA_real_, nrow(df))
        out[as.integer(rownames(df2))] <- stats::prcomp(df2, center = TRUE, scale. = TRUE)$x[, 1]
        return(out)
      }

      eig <- suppressWarnings(WGCNA::moduleEigengenes(df2, cols)$eigengenes)
      if (is.null(eig) || ncol(eig) < 1L) {
        out <- rep(NA_real_, nrow(df))
        out[as.integer(rownames(df2))] <- stats::prcomp(df2, center = TRUE, scale. = TRUE)$x[, 1]
        return(out)
      }

      out <- rep(NA_real_, nrow(df))
      out[as.integer(rownames(df2))] <- as.numeric(eig[, 1])
      return(out)
    }

    stop("Unknown method: ", method)
  }

  network_scalar <- function(g, agg) {
    if (!requireNamespace("igraph", quietly = TRUE)) stop("`graph` requires {igraph} (Suggested).")
    g <- igraph::simplify(g)

    comm <- igraph::cluster_louvain(g)
    Q <- igraph::modularity(comm)

    C_vec <- igraph::transitivity(g, type = "local", isolates = "zero")
    C <- mean(C_vec, na.rm = TRUE)

    D <- igraph::distances(g)
    D[!is.finite(D)] <- NA_real_
    Eglob <- 1 / mean(D, na.rm = TRUE)

    H <- igraph::centralization.degree(g)$centralization

    Qs <- pmin(pmax(Q, 0), 1)
    Cs <- pmin(pmax(C, 0), 1)
    Hs <- pmin(pmax(H, 0), 1)
    Es <- scale_01(Eglob)

    metrics <- c(Qs, Cs, Es, (1 - Hs))
    if (agg == "mean") {
      return(mean(metrics, na.rm = TRUE))
    }
    sqrt(mean(metrics, na.rm = TRUE))
  }

  roll_apply <- function(x, window, fun, align) {
    nloc <- length(x)
    out <- rep(NA_real_, nloc)
    for (i in seq_len(nloc)) {
      if (align == "right") {
        idx <- (i - window + 1):i
      } else if (align == "left") {
        idx <- i:(i + window - 1)
      } else {
        half <- floor(window / 2)
        idx <- (i - half):(i - half + window - 1)
      }
      idx <- idx[idx >= 1 & idx <= nloc]
      if (length(idx) < window) next
      out[i] <- fun(x[idx])
    }
    out
  }

  event_features <- function(z, event, baseline_label, recovery_labels) {
    ok <- is.finite(z) & !is.na(event)
    if (sum(ok) < 3) {
      return(list(baseline = NA_real_, pulse = NA_real_, recovery = NA_real_))
    }

    z <- z[ok]
    e <- event[ok]
    base_idx <- which(e == baseline_label)
    rec_idx <- which(e %in% recovery_labels)
    stress_idx <- which(e != baseline_label)

    if (length(base_idx) < 1) {
      return(list(baseline = NA_real_, pulse = NA_real_, recovery = NA_real_))
    }
    baseline <- mean(z[base_idx])

    pulse <- if (length(stress_idx) >= 1) max(abs(z[stress_idx] - baseline)) else NA_real_
    recovery <- if (length(rec_idx) >= 1) 1 - abs(mean(z[rec_idx]) - baseline) else NA_real_

    list(baseline = baseline, pulse = pulse, recovery = recovery)
  }

  coupling_term <- function(P, S, M, fun) {
    if (fun == "geometric_mean") (P * S * M)^(1 / 3) else (1 - stats::sd(c(P, S, M)))
  }

  # ---- impute (before latent extraction) ----------------------------------
  ROS_flux <- impute_median_by(ROS_flux, id, scale_by)
  Eh_stability <- impute_median_by(Eh_stability, id, scale_by)
  if (!is.null(micro_data)) micro_data <- impute_median_by(micro_data, id, scale_by)

  # ---- latents (raw + oriented) -------------------------------------------
  phys_raw <- orient_latent(latent_dimension(ROS_flux, method_phys), ROS_flux, direction_phys, direction_anchor_phys)
  soil_raw <- orient_latent(latent_dimension(Eh_stability, method_soil), Eh_stability, direction_soil, direction_anchor_soil)

  micro_abund_raw <- rep(NA_real_, n)
  if (!is.null(micro_data)) {
    micro_abund_raw <- orient_latent(latent_dimension(micro_data, method_micro), micro_data, direction_micro, direction_anchor_micro)
  }

  micro_net_raw <- rep(NA_real_, n)
  if (!is.null(graph)) {
    if (inherits(graph, "igraph")) {
      micro_net_raw <- rep(network_scalar(graph, network_agg), n)
    } else if (is.list(graph)) {
      if (length(graph) != n) stop("If `graph` is a list, it must have length nrow(ROS_flux).")
      micro_net_raw <- vapply(graph, function(g) network_scalar(g, network_agg), numeric(1))
    } else {
      stop("`graph` must be NULL, an igraph object, or a list of igraph objects.")
    }
  }

  if (all(is.na(micro_abund_raw)) && all(is.na(micro_net_raw))) {
    stop("Provide microbial input via `micro_data` and/or `graph`.")
  }

  # ---- scale to [0,1] (optional within strata) ----------------------------
  phys <- scale_01_by(phys_raw, id, scale_by)
  soil <- scale_01_by(soil_raw, id, scale_by)
  micro_abund <- scale_01_by(micro_abund_raw, id, scale_by)
  micro_net <- scale_01_by(micro_net_raw, id, scale_by)

  # ---- blend microbial domain ---------------------------------------------
  micro <- if (all(is.na(micro_abund))) {
    micro_net
  } else if (all(is.na(micro_net))) {
    micro_abund
  } else {
    scale_01(alpha_micro * micro_abund + (1 - alpha_micro) * micro_net)
  }

  # ---- snapshot RRI -------------------------------------------------------
  base_rri <- w1 * phys + w2 * soil + w3 * micro

  if (isTRUE(add_coupling) && coupling_weight > 0) {
    ww <- c(w1, w2, w3)
    ww <- ww / sum(ww) * (1 - coupling_weight)

    coup <- mapply(
      coupling_term,
      P = phys,
      S = soil,
      M = micro,
      MoreArgs = list(fun = coupling_fun)
    )

    base_rri <- ww[1] * phys + ww[2] * soil + ww[3] * micro + coupling_weight * coup
  }

  rri <- scale_01(base_rri)

  row_scores <- data.frame(
    Physio = phys,
    Soil = soil,
    Micro = micro,
    RRI = rri,
    Micro_abundance = micro_abund,
    Micro_network = micro_net
  )

  # ---- compositional projection (ternary-ready) ---------------------------
  comp <- row_scores[, c("Physio", "Soil", "Micro")]
  s <- rowSums(comp)
  s[s == 0] <- NA_real_ # undefined composition when all are zero
  comp <- comp / s

  row_scores_comp <- cbind(comp, RRI = row_scores$RRI)
  rri_index <- mean(row_scores$RRI, na.rm = TRUE)
  attr(row_scores_comp, "RRI_index") <- rri_index

  # ---- dynamic modes -------------------------------------------------------

  dyn_scores <- NULL

  if (mode != "snapshot") {
    id2 <- id
    id2[[time_col]] <- suppressWarnings(as.numeric(id2[[time_col]]))

    key <- interaction(id2[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
    ord <- order(key, id2[[time_col]])

    rs_ord <- row_scores[ord, , drop = FALSE]
    id_ord <- id2[ord, , drop = FALSE]
    key_ord <- key[ord]

    if (mode == "rolling") {
      stab_fun <- function(x) 1 - stats::sd(x)

      P_roll <- S_roll <- M_roll <- rep(NA_real_, nrow(rs_ord))
      P_stab <- S_stab <- M_stab <- rep(NA_real_, nrow(rs_ord))

      for (k in unique(key_ord)) {
        idx <- which(key_ord == k)
        P_roll[idx] <- roll_apply(rs_ord$Physio[idx], window, mean, align)
        S_roll[idx] <- roll_apply(rs_ord$Soil[idx], window, mean, align)
        M_roll[idx] <- roll_apply(rs_ord$Micro[idx], window, mean, align)

        P_stab[idx] <- roll_apply(rs_ord$Physio[idx], window, stab_fun, align)
        S_stab[idx] <- roll_apply(rs_ord$Soil[idx], window, stab_fun, align)
        M_stab[idx] <- roll_apply(rs_ord$Micro[idx], window, stab_fun, align)
      }

      P_dyn <- scale_01(0.5 * P_roll + 0.5 * P_stab)
      S_dyn <- scale_01(0.5 * S_roll + 0.5 * S_stab)
      M_dyn <- scale_01(0.5 * M_roll + 0.5 * M_stab)

      rri_dyn <- scale_01(w1 * P_dyn + w2 * S_dyn + w3 * M_dyn)

      dyn_scores <- data.frame(
        P_level = P_roll,
        P_stability = P_stab,
        S_level = S_roll,
        S_stability = S_stab,
        M_level = M_roll,
        M_stability = M_stab,
        Physio_dyn = P_dyn,
        Soil_dyn = S_dyn,
        Micro_dyn = M_dyn,
        RRI_dyn = rri_dyn
      )
    } else if (mode == "event") {
      if (is.null(event_col) || !event_col %in% names(id_ord)) stop("mode = 'event' requires `event_col` in `id`.")
      e <- id_ord[[event_col]]

      out_list <- vector("list", length(unique(key_ord)))
      names(out_list) <- unique(key_ord)

      for (k in unique(key_ord)) {
        idx <- which(key_ord == k)

        P_feat <- event_features(rs_ord$Physio[idx], e[idx], baseline_label, recovery_labels)
        S_feat <- event_features(rs_ord$Soil[idx], e[idx], baseline_label, recovery_labels)
        M_feat <- event_features(rs_ord$Micro[idx], e[idx], baseline_label, recovery_labels)

        P_res <- scale_01(1 - P_feat$pulse) * scale_01(P_feat$recovery)
        S_res <- scale_01(1 - S_feat$pulse) * scale_01(S_feat$recovery)
        M_res <- scale_01(1 - M_feat$pulse) * scale_01(M_feat$recovery)

        rri_e <- scale_01(w1 * P_res + w2 * S_res + w3 * M_res)
        gdat <- id_ord[idx[1], group_cols, drop = FALSE]

        out_list[[k]] <- cbind(
          gdat,
          P_baseline = P_feat$baseline,
          P_pulse = P_feat$pulse,
          P_recovery = P_feat$recovery,
          S_baseline = S_feat$baseline,
          S_pulse = S_feat$pulse,
          S_recovery = S_feat$recovery,
          M_baseline = M_feat$baseline,
          M_pulse = M_feat$pulse,
          M_recovery = M_feat$recovery,
          Physio_event = P_res,
          Soil_event = S_res,
          Micro_event = M_res,
          RRI_event = rri_e
        )
      }

      dyn_scores <- do.call(rbind, out_list)
      rownames(dyn_scores) <- NULL
    }
  }

  out <- list(
    row_scores = row_scores,
    row_scores_comp = row_scores_comp,
    dyn_scores = dyn_scores,
    meta = list(
      mode = mode,
      time_col = time_col,
      group_cols = group_cols,
      scale_by = scale_by,
      method = list(phys = method_phys, soil = method_soil, micro = method_micro),
      direction = list(phys = direction_phys, soil = direction_soil, micro = direction_micro),
      alpha_micro = alpha_micro,
      weights = list(
        w1 = w1, w2 = w2, w3 = w3,
        add_coupling = add_coupling,
        coupling_weight = coupling_weight,
        coupling_fun = coupling_fun
      ),
      rri_index = rri_index
    )
  )

  class(out) <- "RRI"
  out
}
