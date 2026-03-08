# Simulate a Holobiont Redox Dataset with Spatio-Temporal Structure

Generate a synthetic plant–soil–microbiome dataset with (i) a latent
redox state evolving over time, (ii) a mid-season disturbance pulse,
(iii) partial cross-domain decoupling, (iv) configurable microbial
sparsity (zero inflation), and (v) missing-not-at-random (MNAR) Eh
dropout under strongly reduced states.

The simulator is designed for benchmarking, validation, and
demonstration of RedoxRRI workflows. It is not intended to represent any
real ecosystem.

## Usage

``` r
simulate_redox_holobiont(
  n_plot = 4,
  n_depth = 2,
  n_plant = 5,
  n_time = 15,
  p_micro = 60,
  seed = 123,
  include_graph = FALSE,
  depth_labels = NULL,
  disturbance_strength = 0.4,
  disturbance_center = NULL,
  disturbance_width = 0.1,
  seasonal_amp = 0.08,
  seasonal_phase = 0,
  stochastic_reassembly = TRUE,
  decoupling = 0.3,
  zero_inflation = 0.4,
  MNAR_strength = 0.4,
  Eh_dropout_threshold = 0.25,
  micro_mean = 8,
  micro_slope = 3,
  micro_lambda_min = 1e-08,
  micro_lambda_max = 1e+06
)
```

## Arguments

- n_plot:

  Integer \>= 1. Number of plots.

- n_depth:

  Integer \>= 1. Number of soil depth strata.

- n_plant:

  Integer \>= 1. Number of plants per plot–depth.

- n_time:

  Integer \>= 2. Number of time points per plant.

- p_micro:

  Integer \>= 1. Number of microbial features.

- seed:

  Optional integer seed. Use `NULL` for stochastic runs.

- include_graph:

  Logical. If `TRUE` and igraph is available, return a simulated network
  as `graph`.

- depth_labels:

  Optional character vector of length `n_depth`. If `NULL`, depths are
  labeled `"D1", "D2", ...`.

- disturbance_strength:

  Numeric in `[0, 1]`. Magnitude of disturbance pulse.

- disturbance_center:

  Optional numeric. If `NULL`, centered at median time.

- disturbance_width:

  Numeric \> 0. Pulse width as fraction of season length.

- seasonal_amp:

  Numeric \>= 0. Amplitude of seasonal forcing.

- seasonal_phase:

  Numeric. Phase shift for seasonal sine wave (radians).

- stochastic_reassembly:

  Logical. If `TRUE`, microbial counts include a stochastic reassembly
  component (mixture of intensities).

- decoupling:

  Numeric in `[0, 1]`. Cross-domain decoupling (0=tightly coupled).

- zero_inflation:

  Numeric in `[0, 1]`. Probability a microbial entry is set to 0.

- MNAR_strength:

  Numeric in `[0, 1]`. Strength of MNAR Eh dropout in reduced states.

- Eh_dropout_threshold:

  Numeric in `[0, 1]`. Latent threshold below which dropout can occur.

- micro_mean:

  Numeric \> 0. Baseline mean intensity for microbial counts.

- micro_slope:

  Numeric \>= 0. Coupling strength between latent redox state and
  microbial intensity.

- micro_lambda_min:

  Numeric \> 0. Lower bound for microbial Poisson intensity.

- micro_lambda_max:

  Numeric \> 0. Upper bound for microbial Poisson intensity.

## Value

A named list with:

- `id`: Data frame of experimental design (plot, depth, plant_id, time).

- `ROS_flux`: Plant physiological indicators (data frame).

- `Eh_stability`: Soil redox chemistry variables (data frame; Eh may
  include NA).

- `micro_data`: Microbial feature matrix (data frame; non-negative
  integers).

- `latent_truth`: Underlying latent redox state in `[0, 1]`.

- `graph`: Optional `igraph` object if requested and available.

## Examples

``` r
# ---- Simulate small holobiont dataset ----
sim <- simulate_redox_holobiont(
  n_plot = 2,
  n_depth = 1,
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

# ---- Correlation with known latent truth (simulation only) ----
cor_val <- rri_latent_correlation(
  res,
  latent_truth = sim$latent_truth,
  method = "pearson"
)
```
