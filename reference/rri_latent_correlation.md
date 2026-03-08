# Latent Recovery Correlation for RedoxRRI

Computes the Pearson correlation between the per-sample Redox Resilience
Index (RRI) and a known latent truth vector. This function is intended
for simulation-based validation only and should not be interpreted as
predictive accuracy.

## Usage

``` r
rri_latent_correlation(
  res,
  latent_truth,
  method = c("pearson", "spearman", "kendall")
)
```

## Arguments

- res:

  An object returned by
  [`rri_pipeline_st()`](https://mghotbi.github.io/RedoxRRI/reference/rri_pipeline_st.md).

- latent_truth:

  Numeric vector of true latent values.

- method:

  Correlation method. One of `"pearson"`, `"spearman"`, or `"kendall"`.

## Value

A single numeric correlation coefficient.

## Details

This function is designed for simulation benchmarking. In empirical
datasets, no latent truth exists and this metric should not be used.
