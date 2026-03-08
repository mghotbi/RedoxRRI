# Holobiont Redox Resilience Index (RRI) with Spatio-Temporal Dynamics

Computes a holobiont-level Redox Resilience Index (RRI) by integrating
plant physiological traits, soil redox chemistry, and microbial
resilience into a unified, directionally identifiable index. The
framework supports static (snapshot), rolling-window, and event-based
resilience modes, optional compositional geometry (clr), multiblock
reduction via MFA, and covariance-based compensation.

## Usage

``` r
rri_pipeline_st(
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
  reducer = c("per_domain", "mfa"),
  scaling = c("minmax_legacy", "pnorm"),
  comp_space = c("closure_legacy", "clr"),
  ref_stats = NULL,
  add_compensation = FALSE,
  compensation_weight = 0
)
```

## Arguments

- ROS_flux:

  Data frame of plant physiological variables (rows = samples).

- Eh_stability:

  Data frame of soil redox chemistry variables (rows = samples).

- micro_data:

  Optional data frame of microbial abundance or functional features.

- graph:

  Optional `igraph` object or list of `igraph` objects representing
  microbial network structure.

- id:

  Optional data frame describing experimental design (same number of
  rows as inputs).

- time_col:

  Optional character. Name of time column in `id`.

- group_cols:

  Optional character vector of grouping variables in `id`.

- mode:

  Character. One of `"snapshot"`, `"rolling"`, or `"event"`.

- window:

  Integer \>= 2. Rolling window size (for mode = "rolling").

- align:

  Character. Alignment rule for rolling window: `"right"`, `"center"`,
  or `"left"`.

- event_col:

  Optional character. Column in `id` identifying event phases.

- baseline_label:

  Character. Label identifying baseline phase.

- recovery_labels:

  Character vector identifying recovery phases.

- alpha_micro:

  Numeric between 0 and 1 controlling blending of microbial abundance
  and network components.

- method_phys:

  Character. Reduction method for plant block.

- method_soil:

  Character. Reduction method for soil block.

- method_micro:

  Character. Reduction method for microbial block.

- direction_phys:

  Character. Orientation rule for plant latent dimension.

- direction_soil:

  Character. Orientation rule for soil latent dimension.

- direction_micro:

  Character. Orientation rule for microbial latent dimension.

- direction_anchor_phys:

  Optional character. Anchor variable for plant orientation.

- direction_anchor_soil:

  Optional character. Anchor variable for soil orientation.

- direction_anchor_micro:

  Optional character. Anchor variable for microbial orientation.

- scale_by:

  Optional character vector of grouping variables used for scaling.

- network_agg:

  Character. Network aggregation method: `"equation"` or `"mean"`.

- w1:

  Numeric weight for plant domain.

- w2:

  Numeric weight for soil domain.

- w3:

  Numeric weight for microbial domain. Must sum with w1 and w2 to 1.

- add_coupling:

  Logical. If TRUE, adds cross-domain coherence term.

- coupling_weight:

  Numeric between 0 and 1 controlling weight of coupling term.

- coupling_fun:

  Character. Coupling function: `"geometric_mean"` or `"agreement"`.

- norm_method:

  Optional character. If provided, overrides block-specific methods.

- reducer:

  Character. Reduction strategy: `"per_domain"` or `"mfa"`.

- scaling:

  Character. Scaling rule: `"minmax_legacy"` or `"pnorm"`.

- comp_space:

  Character. Compositional projection method: `"closure_legacy"` or
  `"clr"`.

- ref_stats:

  Optional list of reference statistics used for scaling.

- add_compensation:

  Logical. If TRUE, includes covariance-based compensation term.

- compensation_weight:

  Numeric between 0 and 1 controlling compensation weight.

## Value

A list of class `"RRI"` containing:

- `row_scores`: Raw domain and RRI values.

- `row_scores_comp`: Compositional domain scores and RRI.

- `dyn_scores`: Dynamic resilience metrics (if applicable).

- `meta`: Metadata describing model configuration.

## Details

When `reducer = "mfa"`, blocks are integrated using FactoMineR multiple
factor analysis. If partial coordinates are unavailable, the function
safely falls back to per-domain reduction.

When `comp_space = "clr"`, domain scores are projected into Aitchison
geometry using centered log-ratio transformation and returned in simplex
form for ternary visualization.

## Examples

``` r
# ---- Simulate small holobiont dataset ----
sim <- simulate_redox_holobiont(
  n_plot = 10,
  n_depth = 10,
  n_plant = 4,
  n_time = 8,
  p_micro = 20,
  seed = 1234
)

# ---- Snapshot RRI computation ----
res <- rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data,
  id = sim$id,
  reducer = "per_domain",
  scaling = "pnorm"
)

# Per-sample domain scores and RRI
head(res$row_scores)
#>       Physio       Soil     Micro        RRI Micro_abundance Micro_network
#> 1 0.12645864 0.27279357 0.2549676 0.16626299       0.2549676            NA
#> 2 0.06540984 0.07158139 0.2339642 0.09034629       0.2339642            NA
#> 3 0.32618020 0.23589246 0.2764347 0.24136082       0.2764347            NA
#> 4 0.24398316 0.04555534 0.2174682 0.13056096       0.2174682            NA
#> 5 0.09059057 0.08666867 0.1287014 0.08393124       0.1287014            NA
#> 6 0.07421309 0.27507927 0.1888545 0.13473642       0.1888545            NA
#>   Micro_mfa
#> 1        NA
#> 2        NA
#> 3        NA
#> 4        NA
#> 5        NA
#> 6        NA

# Compositional (ternary-ready) allocation
head(res$row_scores_comp)
#>      Physio       Soil     Micro        RRI
#> 1 0.1932969 0.41697538 0.3897278 0.16626299
#> 2 0.1763280 0.19296492 0.6307071 0.09034629
#> 3 0.3890010 0.28132425 0.3296748 0.24136082
#> 4 0.4812228 0.08985155 0.4289257 0.13056096
#> 5 0.2960857 0.28326743 0.4206469 0.08393124
#> 6 0.1379049 0.51116024 0.3509349 0.13473642

# ---- Rolling dynamic mode example ----
res_roll <- rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data,
  id = sim$id,
  mode = "rolling",
  time_col = "time",
  group_cols = c("plot", "depth", "plant_id"),
  window = 2
)

# Note: The first (window - 1) rows per group are NA
# due to right-aligned rolling windows.
head(res_roll$dyn_scores)
#>     P_level P_stability   S_level S_stability    M_level M_stability Physio_dyn
#> 1        NA          NA        NA          NA         NA          NA         NA
#> 2 0.1132156   0.8398890 0.1664783   0.8500662 0.08284815   0.9880108  0.1327787
#> 3 0.2129241   0.9808981 0.1399550   0.8875758 0.09741222   0.9913925  0.3322787
#> 4 0.4821401   0.6001690 0.5150332   0.5819835 0.55174934   0.3660778  0.2398599
#> 5 0.8824317   0.8337333 0.9053078   0.8660850 0.78810327   0.7003328  0.7651818
#> 6 0.6677582   0.5301391 0.8094827   0.7305679 0.59126873   0.9786988  0.3356560
#>    Soil_dyn Micro_dyn   RRI_dyn
#> 1        NA        NA        NA
#> 2 0.1853902 0.2306539 0.1769542
#> 3 0.1944970 0.2455451 0.2676863
#> 4 0.2520957 0.1036692 0.2129845
#> 5 0.8111031 0.5771562 0.7614496
#> 6 0.6193377 0.6448105 0.5291372

# System-level mean RRI
attr(res$row_scores_comp, "RRI_index")
#> [1] 0.4792014
```
