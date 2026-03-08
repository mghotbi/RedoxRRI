# Cross-Domain Compensation Index

Computes a scale-invariant cross-domain compensation index based on mean
pairwise correlations among Physio, Soil, and Micro domain scores.

## Usage

``` r
rri_compensation_index(res)
```

## Arguments

- res:

  An object returned by
  [`rri_pipeline_st()`](https://mghotbi.github.io/RedoxRRI/reference/rri_pipeline_st.md).

## Value

A single numeric compensation index in \\\[-1, 1\]\\. Positive values
indicate compensatory dynamics (negative cross-domain correlations),
while negative values indicate synchrony.

## Details

The compensation index is defined as:

\$\$ \mathrm{Comp} = -\frac{2}{k(k-1)} \sum\_{i\<j} \mathrm{Cor}(D_i,
D_j) \$\$

where \\k = 3\\ and \\D_i\\ are domain-level scores.

Using correlations rather than covariances ensures scale invariance and
boundedness across datasets.

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

# ---- Compute cross-domain compensation index ----
comp <- rri_compensation_index(res)
```
