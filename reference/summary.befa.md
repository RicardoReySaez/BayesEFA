# Summary and Print Methods for BEFA Objects

Methods to summarize and display results from Bayesian Exploratory
Factor Analysis.

Computes posterior summaries for a fitted Bayesian EFA model, including
point estimates, credible intervals, and convergence diagnostics for
factor loadings.

## Usage

``` r
# S3 method for class 'befa'
summary(
  object,
  probs = c(0.025, 0.975),
  cutoff = 0,
  sort = FALSE,
  signif_stars = FALSE,
  ...
)
```

## Arguments

- object:

  A `befa` object returned by
  [`befa()`](https://ricardoreysaez.github.io/BayesEFA/reference/befa.md).

- probs:

  Numeric vector of length 2. Quantiles for credible intervals (default:
  0.025, 0.975).

- cutoff:

  Numeric. Loadings with absolute value below this threshold are hidden
  in print output.

- sort:

  Logical. If TRUE, items are sorted by their primary factor loading.

- signif_stars:

  Logical. If TRUE, marks loadings whose CI excludes zero with an
  asterisk.

- ...:

  Ignored.

## Value

A `summary.befa` object containing:

- `estimates`: Data frame with loadings, h2, u2, and diagnostics per
  item.

- `sig_matrix`: Logical matrix indicating which loadings have CIs
  excluding zero.

- `probs`: The quantile probabilities used for credible intervals.

- `phi`: Factor correlation matrix (identity for orthogonal rotations).

- `variance`: List with per-factor and total explained variance
  proportions.

- `fit_indices`: Fit measures from the `befa` object (if available).

- `reliability`: Reliability estimates from the `befa` object (if
  available).

- `header`: List with model metadata (model_type, n_factors, N,
  rotation).

- `print_options`: List of display options (cutoff, sort, signif_stars).

- `tables`: Pre-built display tables for loadings, fit measures, and
  reliability.

## Details

The summary extracts posterior draws for loadings (Lambda), computes
means and quantiles, and reports MCMC diagnostics (Rhat, ESS).
Communalities (h2) and uniquenesses (u2) are derived from the squared
loadings.

## Examples

``` r
if (FALSE) { # \dontrun{
befa_fit <- befa(
  data = HS_data, n_factors = 3, model = "cor",
  iter = 500, chains = 2, seed = 123
)

# Get and print the summary
befa_summary <- summary(befa_fit, sort = TRUE, signif_stars = TRUE)
print(befa_summary)
} # }
```
