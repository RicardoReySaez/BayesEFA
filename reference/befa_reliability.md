# Bayesian Reliability Estimation (McDonald's Omega)

Computes model-based reliability coefficients (McDonald's Omega) from
the posterior distribution of a fitted BEFA model. By evaluating the
reliability analytically at each MCMC draw, it provides full posterior
distributions (and credible intervals) for both total-score and subscale
reliabilities, avoiding the need for asymptotic approximations.

## Usage

``` r
befa_reliability(object, probs = c(0.025, 0.975))
```

## Arguments

- object:

  A `befa` object returned by
  [`befa()`](https://ricardoreysaez.github.io/BayesEFA/reference/befa.md).

- probs:

  Numeric vector of length 2. Quantiles used to compute the credible
  intervals (default: `c(0.025, 0.975)` for a 95% interval).

## Value

A list of class `befa_reliability` containing:

- `omega_draws`: A `draws_array` (iterations \\\times\\ chains
  \\\times\\ P) containing the posterior draws of all reliability
  coefficients (`omega_total`, `omega_F1`, ..., `omega_FM`).

- `summary`: A `draws_summary` data frame computed via
  [`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html).

## Details

### Why McDonald's Omega?

Unlike Cronbach's \\\alpha\\, which strictly assumes tau-equivalence
(equal factor loadings across all items), McDonald's \\\omega\\
correctly accounts for the congeneric factor structure estimated in
Exploratory Factor Analysis (EFA). Furthermore, by computing the
coefficients draw-by-draw from the posterior distribution, this function
naturally propagates parameter uncertainty, yielding full credible
intervals for reliability estimates.

### Mathematical Formulation

Let \\J\\ be the number of items, \\M\\ the number of factors,
\\\mathbf{\Lambda}\\ the \\J \times M\\ loading matrix,
\\\mathbf{\Psi}\\ the diagonal matrix of unique variances
(uniquenesses), and \\\mathbf{\Sigma}\\ the model-implied covariance
matrix. \\\mathbf{1}\\ is a vector of ones of length \\J\\.

- **Omega Total (\\\omega_t\\)**: Represents the proportion of variance
  in the unit-weighted total score (sum of all items) attributable to
  *all* common factors. \$\$\omega_t = \frac{\sum\_{m=1}^M
  (\sum\_{j=1}^J \lambda\_{jm})^2}{\mathbf{1}^\top \cdot \mathbf{\Sigma}
  \cdot \mathbf{1}}\$\$

- **Omega Subscale (\\\omega\_{s,m}\\)**: Represents the reliability of
  the sum score when isolating the variance attributable to a specific
  factor \\m\\, treating the unique variances as measurement error.
  \$\$\omega\_{s,m} = \frac{(\sum\_{j=1}^J
  \lambda\_{jm})^2}{(\sum\_{j=1}^J \lambda\_{jm})^2 + \sum\_{j=1}^J
  \psi_j}\$\$

## References

McDonald, R. P. (1999). *Test Theory: A Unified Treatment*. Lawrence
Erlbaum Associates.

Zinbarg, R. E., Revelle, W., Yovel, I., & Li, W. (2005). Cronbach's
\\\alpha\\, Revelle's \\\beta\\, and McDonald's \\\omega_H\\: Their
relations with each other and two alternative conceptualizations of
reliability. *Psychometrika, 70*(1), 123-133.
<https://doi.org/10.1007/s11336-003-0974-7>

## Examples

``` r
if (FALSE) { # \dontrun{
# ------------------------------------------------------------------------ #
#    1. Bayesian Reliability estimates after fitting Bayesian EFA model    #
# ------------------------------------------------------------------------ #

# Fit Bayesian EFA model
befa_fit <- befa(
  data = HS_data,
  n_factors = 3,
  factor_scores = FALSE,
  compute_fit_indices = FALSE,
  compute_reliability = FALSE,
  backend = "rstan",
  seed = 17,
  chains = 4,
  parallel_chains = 4
)

# Compute Bayesian Reliability estimates from a fitted befa object
bayesian_omega <- befa_reliability(befa_fit, probs = c(.025, .975))

# Posterior summaries
bayesian_omega$summary

# # A tibble: 4 × 10
#   variable     mean median     sd    mad    q5   q95  rhat ess_bulk ess_tail
#   <chr>       <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
# 1 omega_total 0.843  0.844 0.0106 0.0105 0.825 0.860 1.000    4659.    3628.
# 2 omega_F1    0.598  0.600 0.0380 0.0377 0.533 0.657 1.00     4057.    3435.
# 3 omega_F2    0.727  0.728 0.0185 0.0181 0.697 0.755 1.000    4429.    3552.
# 4 omega_F3    0.545  0.550 0.0442 0.0420 0.463 0.609 1.00     3489.    3584.

# Posterior draws (transform to matrix)
posterior::as_draws_matrix(bayesian_omega$omega_draws)

# # A draws_matrix: 1000 iterations, 4 chains, and 4 variables
#     variable
# draw omega_total omega_F1 omega_F2 omega_F3
#   1         0.85     0.60     0.73     0.59
#   2         0.85     0.60     0.73     0.56
#   3         0.85     0.61     0.75     0.56
#   4         0.84     0.55     0.73     0.53
#   5         0.86     0.66     0.74     0.55
#   6         0.84     0.58     0.70     0.58
#   7         0.85     0.64     0.73     0.51
#   8         0.86     0.69     0.74     0.46
#   9         0.87     0.56     0.78     0.63
#   10        0.83     0.61     0.72     0.46
# # ... with 3990 more draws

# We can plot the posterior distribution of reliability estimates
par(mfrow = c(2, 2))
hist(posterior::as_draws_matrix(bayesian_omega$omega_draws)[, "omega_total"],
     breaks = 100, col = "grey90",
     main = "Histogram of omega index: Full scale",
     xlab = "Omega"
)
abline(
  v = bayesian_omega$summary$mean[1], col = "firebrick2",
  lwd = 3, lty = 2
)
hist(posterior::as_draws_matrix(bayesian_omega$omega_draws)[, "omega_F1"],
     breaks = 100, col = "grey90",
     main = "Histogram of omega index: Factor 1",
     xlab = "Omega"
)
abline(
  v = bayesian_omega$summary$mean[2], col = "firebrick2",
  lwd = 3, lty = 2
)
hist(posterior::as_draws_matrix(bayesian_omega$omega_draws)[, "omega_F2"],
     breaks = 100, col = "grey90",
     main = "Histogram of omega index: Factor 2",
     xlab = "Omega"
)
abline(
  v = bayesian_omega$summary$mean[3], col = "firebrick2",
  lwd = 3, lty = 2
)
hist(posterior::as_draws_matrix(bayesian_omega$omega_draws)[, "omega_F3"],
     breaks = 100, col = "grey90",
     main = "Histogram of omega index: Factor 3",
     xlab = "Omega"
)
abline(
  v = bayesian_omega$summary$mean[4], col = "firebrick2",
  lwd = 3, lty = 2
)
par(mfrow = c(1, 1))

# -----------------------------------------------------------------------------

# ----------------------------------------------------------------------- #
#    2. Bayesian Reliability estimates when fitting Bayesian EFA model    #
# ----------------------------------------------------------------------- #

# Fit Bayesian model with compute_fit_indices = TRUE
befa_fit <- befa(
  data = HS_data,
  n_factors = 3,
  factor_scores = FALSE,
  compute_fit_indices = FALSE,
  compute_reliability = TRUE,
  backend = "rstan",
  seed = 17,
  chains = 4,
  parallel_chains = 4
)

# befa_reliability output is inside the "reliability" object
# befa_fit$reliability

# -----------------------------------------------------------------------------
} # }
```
