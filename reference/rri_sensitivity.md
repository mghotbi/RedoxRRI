# Sensitivity Analysis for RedoxRRI

Evaluates robustness of RRI rankings under alternative domain
aggregation weights using already computed domain scores.

## Usage

``` r
rri_sensitivity(res, weight_grid = seq(0.2, 0.6, by = 0.1))
```

## Arguments

- res:

  An object returned by
  [`rri_pipeline_st()`](https://mghotbi.github.io/RedoxRRI/reference/rri_pipeline_st.md).

- weight_grid:

  Numeric vector of Physio weights in (0,1).

## Value

A data frame containing tested weights and Spearman rank correlations
relative to the baseline RRI.

## Details

This function perturbs only the aggregation step and does not recompute
latent domain scores. It therefore evaluates stability of domain
integration rather than reducer sensitivity.

## Examples

``` r
# ---- Simulate a small holobiont dataset ----
sim <- simulate_redox_holobiont(
  n_plot = 2,
  n_depth = 4,
  n_plant = 2,
  n_time = 8,
  p_micro = 20,
  seed = 1
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

# ---- Evaluate aggregation sensitivity ----
sens <- rri_sensitivity(
  res,
  weight_grid = seq(0.3, 0.5, by = 0.1)
)

```
