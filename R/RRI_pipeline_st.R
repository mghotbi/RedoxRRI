#' Holobiont Redox Resilience Index (RRI) with spatio-temporal dynamics
#'
#' Computes a holobiont-level Redox Resilience Index (RRI) by integrating
#' plant physiology, soil redox chemistry, and microbial resilience into a unified,
#' directionally identifiable index with explicit support for spatial and temporal designs.
#'
#' @importFrom stats prcomp lm sd coef cor median cov complete.cases pnorm
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
    norm_method = NULL,
    
    # ---- extensions (backward compatible) ----
    reducer = c("per_domain", "mfa"),
    scaling = c("minmax_legacy", "pnorm"),
    comp_space = c("closure_legacy", "clr"),
    ref_stats = NULL,
    add_compensation = FALSE,
    compensation_weight = 0
) {
  mode <- match.arg(mode)
  align <- match.arg(align)
  network_agg <- match.arg(network_agg)
  coupling_fun <- match.arg(coupling_fun)
  direction_phys <- match.arg(direction_phys)
  direction_soil <- match.arg(direction_soil)
  direction_micro <- match.arg(direction_micro)
  reducer <- match.arg(reducer)
  scaling <- match.arg(scaling)
  comp_space <- match.arg(comp_space)
  
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
  if (!is.finite(compensation_weight) || compensation_weight < 0 || compensation_weight > 1) {
    stop("compensation_weight must be in [0,1].")
  }
  
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
    if (all(is.na(v))) return(rep(NA_real_, length(v)))
    r <- range(v, na.rm = TRUE)
    if (!all(is.finite(r)) || diff(r) == 0) return(rep(0.5, length(v)))
    (v - r[1]) / diff(r)
  }
  
  scale_pnorm <- function(v, stats = NULL) {
    if (all(is.na(v))) return(rep(NA_real_, length(v)))
    mu <- if (is.list(stats) && !is.null(stats$mean)) stats$mean else mean(v, na.rm = TRUE)
    si <- if (is.list(stats) && !is.null(stats$sd)) stats$sd else stats::sd(v, na.rm = TRUE)
    if (!is.finite(si) || si == 0) si <- 1
    stats::pnorm((v - mu) / si)
  }
  
  scale_vec <- function(v, stats = NULL) {
    if (scaling == "minmax_legacy") scale_01(v) else scale_pnorm(v, stats)
  }
  
  scale_vec_by <- function(v, id, by, stats = NULL) {
    if (is.null(by) || is.null(id) || length(by) == 0) {
      return(scale_vec(v, stats))
    }
    key <- interaction(id[, by, drop = FALSE], drop = TRUE, lex.order = TRUE)
    out <- v
    for (k in unique(key)) {
      idx <- which(key == k)
      out[idx] <- scale_vec(v[idx], stats)
    }
    out
  }
  
  orient_latent <- function(latent, df, direction, anchor) {
    if (direction == "higher_is_better") return(latent)
    if (direction == "lower_is_better") return(-latent)
    
    if (is.null(anchor) || !anchor %in% names(df)) return(latent)
    
    a <- df[[anchor]]
    ok <- is.finite(latent) & is.finite(a)
    if (!any(ok)) return(latent)
    
    if (stats::cor(latent[ok], a[ok]) < 0) -latent else latent
  }
  
  latent_dimension <- function(df, method) {
    df <- as_numeric_df(df)
    if (ncol(df) == 1L) return(df[[1L]])
    
    if (method == "mean") return(rowMeans(df, na.rm = TRUE))
    if (method == "scale") return(rowMeans(scale(df), na.rm = TRUE))
    if (method == "pca") return(stats::prcomp(df, center = TRUE, scale. = TRUE)$x[, 1])
    
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
    if (agg == "mean") return(mean(metrics, na.rm = TRUE))
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
    if (sum(ok) < 3) return(list(baseline = NA_real_, pulse = NA_real_, recovery = NA_real_))
    
    z <- z[ok]
    e <- event[ok]
    base_idx <- which(e == baseline_label)
    rec_idx <- which(e %in% recovery_labels)
    stress_idx <- which(e != baseline_label)
    
    if (length(base_idx) < 1) return(list(baseline = NA_real_, pulse = NA_real_, recovery = NA_real_))
    baseline <- mean(z[base_idx])
    
    pulse <- if (length(stress_idx) >= 1) max(abs(z[stress_idx] - baseline)) else NA_real_
    recovery <- if (length(rec_idx) >= 1) 1 - abs(mean(z[rec_idx]) - baseline) else NA_real_
    
    list(baseline = baseline, pulse = pulse, recovery = recovery)
  }
  
  coupling_term <- function(P, S, M, fun) {
    if (fun == "geometric_mean") (P * S * M)^(1 / 3) else (1 - stats::sd(c(P, S, M)))
  }
  
  compensation_index <- function(P, S, M) {
    X <- cbind(P, S, M)
    if (nrow(X) < 3) return(NA_real_)
    C <- stats::cov(X, use = "pairwise.complete.obs")
    if (any(!is.finite(C))) return(NA_real_)
    -sum(C[upper.tri(C)])
  }
  
  clr_transform <- function(X, eps = 1e-8) {
    X <- pmax(as.matrix(X), eps)
    gm <- exp(rowMeans(log(X)))
    log(X / gm)
  }
  clr_to_simplex <- function(CLR) {
    X <- exp(CLR)
    X / rowSums(X)
  }
  
  # ---- MFA partial extractor (ONLY if available; otherwise fallback) -------
  extract_mfa_group_dim1 <- function(mfa, group_names) {
    # Some FactoMineR versions expose partial coords, others don't.
    part <- NULL
    
    if (!is.null(mfa$partial) &&
        !is.null(mfa$partial$ind) &&
        !is.null(mfa$partial$ind$coord)) {
      part <- mfa$partial$ind$coord
    } else if (!is.null(mfa$ind) &&
               !is.null(mfa$ind$coord.partiel)) {
      part <- mfa$ind$coord.partiel
    } else if (!is.null(mfa$ind) &&
               !is.null(mfa$ind$coord.partial)) {
      part <- mfa$ind$coord.partial
    }
    
    if (is.null(part)) return(NULL)
    
    # Case 1: true 3D array [n x ncp x g]
    if (length(dim(part)) == 3) {
      out <- lapply(seq_along(group_names), function(k) as.numeric(part[, 1, k]))
      names(out) <- group_names
      return(out)
    }
    
    # Case 2: data.frame/matrix with per-group Dim1 columns
    part_df <- as.data.frame(part)
    cn <- colnames(part_df)
    
    pick <- function(g) {
      # common FactoMineR naming variants
      pats <- c(
        paste0("^", g, "\\.Dim\\.?1$"),
        paste0("^", g, "\\.Dim\\.1$"),
        paste0("^", g, ".*Dim\\.?1$"),
        paste0("^", g, ".*1$")
      )
      for (p in pats) {
        hit <- which(grepl(p, cn))
        if (length(hit) == 1) return(hit)
      }
      integer(0)
    }
    
    idxs <- lapply(group_names, pick)
    if (all(vapply(idxs, length, integer(1)) == 1L)) {
      out <- lapply(group_names, function(g) as.numeric(part_df[[idxs[[g]]]]))
      names(out) <- group_names
      return(out)
    }
    
    # If it’s only ncp columns (like you saw: ncol=5), it’s NOT partials.
    NULL
  }
  
  # ---- impute (before latent extraction) ----------------------------------
  ROS_flux <- impute_median_by(ROS_flux, id, scale_by)
  Eh_stability <- impute_median_by(Eh_stability, id, scale_by)
  if (!is.null(micro_data)) micro_data <- impute_median_by(micro_data, id, scale_by)
  
  # ---- microbial inputs: abundance + network -------------------------------
  micro_abund_raw <- rep(NA_real_, n)
  if (!is.null(micro_data)) {
    micro_abund_raw <- orient_latent(
      latent_dimension(micro_data, method_micro),
      micro_data,
      direction_micro,
      direction_anchor_micro
    )
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
  
  # ---- Physio & Soil raw scores + MFA (robust) -----------------------------
  mfa_global_dim1 <- rep(NA_real_, n)
  mfa_used_partials <- FALSE
  mfa_fallback_to_per_domain <- FALSE
  
  if (reducer == "mfa") {
    if (!requireNamespace("FactoMineR", quietly = TRUE)) {
      stop("reducer = 'mfa' requires {FactoMineR}.")
    }
    
    micro_block <- NULL
    if (!is.null(micro_data)) micro_block <- as_numeric_df(micro_data)
    if (!all(is.na(micro_net_raw))) {
      micro_net_df <- data.frame(net_resilience = as.numeric(micro_net_raw))
      micro_block <- if (is.null(micro_block)) micro_net_df else cbind(micro_block, micro_net_df)
    }
    if (is.null(micro_block)) stop("MFA requires microbial input via `micro_data` and/or `graph`.")
    
    blocks <- list(
      Physio = as_numeric_df(ROS_flux),
      Soil   = as_numeric_df(Eh_stability),
      Micro  = as_numeric_df(micro_block)
    )
    
    X <- do.call(cbind, blocks)
    groups <- vapply(blocks, ncol, integer(1))
    group_names <- names(groups)
    
    mfa <- FactoMineR::MFA(
      X,
      group = groups,
      type = rep("s", length(groups)),
      name.group = group_names,
      ncp = 5,
      graph = FALSE
    )
    
    if (!is.null(mfa$ind) && !is.null(mfa$ind$coord) && ncol(mfa$ind$coord) >= 1) {
      mfa_global_dim1 <- as.numeric(mfa$ind$coord[, 1])
    }
    
    parts <- extract_mfa_group_dim1(mfa, group_names = group_names)
    
    if (!is.null(parts)) {
      # TRUE partials exist -> safe, non-degenerate ternary
      mfa_used_partials <- TRUE
      phys_raw <- as.numeric(parts$Physio)
      soil_raw <- as.numeric(parts$Soil)
      micro_mfa_raw <- as.numeric(parts$Micro)
      
      phys_raw <- orient_latent(phys_raw, ROS_flux, direction_phys, direction_anchor_phys)
      soil_raw <- orient_latent(soil_raw, Eh_stability, direction_soil, direction_anchor_soil)
      
      micro_anchor_df <- if (!is.null(micro_data)) micro_data else micro_block
      micro_mfa_raw <- orient_latent(micro_mfa_raw, micro_anchor_df, direction_micro, direction_anchor_micro)
      
    } else {
      # No partials available -> fallback to per_domain for domain scores
      mfa_fallback_to_per_domain <- TRUE
      
      phys_raw <- orient_latent(
        latent_dimension(ROS_flux, method_phys),
        ROS_flux,
        direction_phys,
        direction_anchor_phys
      )
      soil_raw <- orient_latent(
        latent_dimension(Eh_stability, method_soil),
        Eh_stability,
        direction_soil,
        direction_anchor_soil
      )
      micro_mfa_raw <- rep(NA_real_, n)
    }
    
  } else {
    phys_raw <- orient_latent(
      latent_dimension(ROS_flux, method_phys),
      ROS_flux,
      direction_phys,
      direction_anchor_phys
    )
    soil_raw <- orient_latent(
      latent_dimension(Eh_stability, method_soil),
      Eh_stability,
      direction_soil,
      direction_anchor_soil
    )
    micro_mfa_raw <- rep(NA_real_, n)
  }
  
  # ---- scale to [0,1] ------------------------------------------------------
  phys <- scale_vec_by(phys_raw, id, scale_by, stats = if (is.list(ref_stats)) ref_stats$phys else NULL)
  soil <- scale_vec_by(soil_raw, id, scale_by, stats = if (is.list(ref_stats)) ref_stats$soil else NULL)
  
  micro_abund <- scale_vec_by(micro_abund_raw, id, scale_by, stats = if (is.list(ref_stats)) ref_stats$micro_abund else NULL)
  micro_net <- scale_vec_by(micro_net_raw, id, scale_by, stats = if (is.list(ref_stats)) ref_stats$micro_net else NULL)
  micro_mfa <- scale_vec_by(micro_mfa_raw, id, scale_by, stats = if (is.list(ref_stats)) ref_stats$micro_mfa else NULL)
  
  # ---- blend microbial domain ---------------------------------------------
  micro <- if (all(is.na(micro_abund))) {
    micro_net
  } else if (all(is.na(micro_net))) {
    micro_abund
  } else {
    scale_vec(alpha_micro * micro_abund + (1 - alpha_micro) * micro_net,
              stats = if (is.list(ref_stats)) ref_stats$micro else NULL
    )
  }
  
  # If MFA partials exist, optionally blend MFA-micro with legacy micro estimate
  if (reducer == "mfa" && isTRUE(mfa_used_partials) && !all(is.na(micro_mfa))) {
    micro <- scale_vec(0.5 * micro + 0.5 * micro_mfa,
                       stats = if (is.list(ref_stats)) ref_stats$micro else NULL
    )
  }
  
  # ---- snapshot RRI -------------------------------------------------------
  base_rri <- w1 * phys + w2 * soil + w3 * micro
  
  if (isTRUE(add_compensation) && compensation_weight > 0) {
    comp_term <- rep(NA_real_, n)
    
    if (!is.null(id) && !is.null(group_cols) && length(group_cols) > 0 && all(group_cols %in% names(id))) {
      keyc <- interaction(id[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
      for (k in unique(keyc)) {
        idx <- which(keyc == k)
        comp_term[idx] <- compensation_index(phys[idx], soil[idx], micro[idx])
      }
    } else {
      comp_term[] <- compensation_index(phys, soil, micro)
    }
    comp_term01 <- scale_vec(
      comp_term,
      stats = if (is.list(ref_stats)) ref_stats$comp else NULL
    )
    
    # Guard against numerical pathologies
    comp_term01[!is.finite(comp_term01)] <- 0
    ww <- c(w1, w2, w3)
    ww <- ww / sum(ww) * (1 - compensation_weight)
    base_rri <- ww[1] * phys + ww[2] * soil + ww[3] * micro + compensation_weight * comp_term01
  }
  
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
  
  rri <- scale_vec(base_rri, stats = if (is.list(ref_stats)) ref_stats$rri else NULL)
  
  row_scores <- data.frame(
    Physio = phys,
    Soil = soil,
    Micro = micro,
    RRI = rri,
    Micro_abundance = micro_abund,
    Micro_network = micro_net,
    Micro_mfa = micro_mfa
  )
  
  # ---- compositional projection (ternary-ready) ---------------------------
  comp <- row_scores[, c("Physio", "Soil", "Micro")]
  
  if (comp_space == "clr") {
    clr <- clr_transform(comp)
    comp_simplex <- clr_to_simplex(clr)
    comp <- as.data.frame(comp_simplex)
    attr(comp, "clr") <- clr
  } else {
    s <- rowSums(comp)
    s[s == 0] <- NA_real_
    comp <- comp / s
  }
  
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
      
      P_dyn <- scale_vec(0.5 * P_roll + 0.5 * P_stab)
      S_dyn <- scale_vec(0.5 * S_roll + 0.5 * S_stab)
      M_dyn <- scale_vec(0.5 * M_roll + 0.5 * M_stab)
      
      rri_dyn <- scale_vec(w1 * P_dyn + w2 * S_dyn + w3 * M_dyn)
      
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
        
        P_res <- scale_vec(1 - P_feat$pulse) * scale_vec(P_feat$recovery)
        S_res <- scale_vec(1 - S_feat$pulse) * scale_vec(S_feat$recovery)
        M_res <- scale_vec(1 - M_feat$pulse) * scale_vec(M_feat$recovery)
        
        rri_e <- scale_vec(w1 * P_res + w2 * S_res + w3 * M_res)
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
      reducer = reducer,
      scaling = scaling,
      comp_space = comp_space,
      add_compensation = add_compensation,
      compensation_weight = compensation_weight,
      weights = list(
        w1 = w1, w2 = w2, w3 = w3,
        add_coupling = add_coupling,
        coupling_weight = coupling_weight,
        coupling_fun = coupling_fun
      ),
      rri_index = rri_index,
      # new but harmless fields:
      mfa_global_dim1 = mfa_global_dim1,
      mfa_used_partials = mfa_used_partials,
      mfa_fallback_to_per_domain = mfa_fallback_to_per_domain
    )
  )
  
  class(out) <- "RRI"
  out
}