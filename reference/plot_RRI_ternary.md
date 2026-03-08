# Ternary Plot of Holobiont Redox Resilience Allocation

Creates a ternary diagram visualising the compositional allocation of
holobiont redox buffering capacity across plant physiology (`Physio`),
soil redox chemistry (`Soil`), and microbial resilience (`Micro`).
Points are filled according to per-sample `RRI` values.

## Usage

``` r
plot_RRI_ternary(
  ternary_df,
  point_size = 5,
  point_alpha = 0.9,
  palette = "plasma",
  show_subtitle = TRUE,
  show_centroid = TRUE,
  centroid_shape = 23,
  centroid_size = 1.4,
  tolerance = 1e-06,
  renormalize = FALSE,
  centroid_method = c("auto", "simplex_mean", "aitchison_mean")
)
```

## Arguments

- ternary_df:

  A data frame containing compositional columns `Physio`, `Soil`,
  `Micro`, and `RRI`.

- point_size:

  Numeric; size of ternary points.

- point_alpha:

  Numeric between 0 and 1 controlling point transparency.

- palette:

  Character; viridis palette option.

- show_subtitle:

  Logical; display system-level RRI mean in subtitle.

- show_centroid:

  Logical; add compositional centroid marker.

- centroid_shape:

  Numeric; ggplot2 shape for centroid marker.

- centroid_size:

  Numeric multiplier for centroid size.

- tolerance:

  Numeric; tolerance used for compositional closure checks.

- renormalize:

  Logical; if TRUE, renormalises rows that do not sum to 1.

- centroid_method:

  Character; one of `"auto"`, `"simplex_mean"`, or `"aitchison_mean"`.

## Value

A `ggtern` object.

## Details

If clr-transformed coordinates are attached as an attribute (`"clr"`),
the centroid can be computed using the Aitchison mean, ensuring
geometric coherence in compositional space. Otherwise, a simplex
arithmetic mean is used.

## Examples

``` r
# ---- Simulate holobiont dataset (snapshot) ----
sim <- simulate_redox_holobiont(
  n_plot = 20,
  n_depth = 30,
  n_plant = 8,
  n_time = 10,
  p_micro = 30,
  seed = 1234
)

# ---- Compute RedoxRRI ----
res <- rri_pipeline_st(
  ROS_flux = sim$ROS_flux,
  Eh_stability = sim$Eh_stability,
  micro_data = sim$micro_data,
  id = sim$id,
  reducer = "per_domain",
  scaling = "pnorm"
)

# ---- Extract compositional scores ----
ternary_df <- res$row_scores_comp

# ---- Plot ternary allocation ----
p <- plot_RRI_ternary(
  ternary_df,
  point_size = 3,
  show_centroid = TRUE
)
#> Warning: Ignoring unknown aesthetics: z
```
