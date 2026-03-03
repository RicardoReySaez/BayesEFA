# Bayesian Fit Indices for Exploratory Factor Analysis

Computes Bayesian analogs of common Structural Equation Modeling (SEM)
fit indices. This implementation follows the **DevM framework** proposed
by Garnier-Villarreal and Jorgensen (2020), which translates traditional
maximum-likelihood fit indices into the Bayesian framework. By
evaluating these indices across the MCMC draws, the function yields full
posterior distributions, allowing for the quantification of sampling
uncertainty (e.g., 95% credible intervals) rather than relying on a
single point estimate.

## Usage

``` r
befa_fit_measures(object, ...)
```

## Arguments

- object:

  A `befa` object returned by
  [`befa()`](https://ricardoreysaez.github.io/BayesEFA/reference/befa.md).

- ...:

  Additional arguments passed to Stan when automatically fitting the
  null model (e.g., `iter`, `chains`, `warmup`, `loo_args`).

## Value

A list of class `befa_fitmeasures` containing:

- `fit_draws`: A `draws_array` (iterations \\\times\\ chains
  \\\times\\ P) containing the posterior draws of all fit indices
  (`Chi2`, `Chi2_Null`, `BRMSEA`, `BGamma`, `Adj_BGamma`, `BMc`, `SRMR`,
  `BCFI`, `BTLI`).

- `summary`: A `draws_summary` data frame computed via
  [`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html).

- `posterior_fit`: A data frame containing the full MCMC posterior draws
  for all \\\chi^2\\-based indices, enabling custom plotting and density
  estimation.

- `loo_object`: The full `loo` object returned by
  [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html).

- `details`: A list with internal model constants (\\p^\*\\, \\pD\\,
  \\N\\, and `chi2_ppp`).

## Details

**The DevM Framework (Deviance at the Posterior Mean)**

In frequentist SEM, fit indices are derived from the exact fit test
statistic (\\\chi^2\\) and the model degrees of freedom (\\df\\). In the
Bayesian DevM framework (Garnier-Villarreal & Jorgensen, 2020), these
quantities are replaced by their Bayesian analogs:

- **Observed Deviance (\\D_i^{obs}\\)**: The model deviance evaluated at
  each posterior draw \\i\\, serving as the Bayesian analog of
  \\\chi^2\\.

- **Effective Number of Parameters (\\pD\\)**: A penalty for model
  complexity. Following the authors' recommendations, this function
  estimates \\pD\\ using the Leave-One-Out Information Criterion
  (\\p\_{loo}\\) via Pareto Smoothed Importance Sampling Leave-One-Out
  (PSIS-LOO, Vehtari et al., 2017).

- **Saturated Moments (\\p^\*\\)**: The number of nonredundant sample
  moments (e.g., \\J \cdot (J + 3) / 2\\ for raw data, or \\J \cdot
  (J + 1) / 2\\ for covariance matrices).

- **Bayesian Degrees of Freedom (\\df_B\\)**: Defined as \\df_B = p^\* -
  pD\\.

- **Bayesian Noncentrality Parameter (\\\lambda_i\\)**: Defined at each
  iteration as \\\lambda_i = D_i^{obs} - p^\*\\.

**Absolute Fit Indices**

Absolute indices evaluate how well the hypothesized model reproduces the
observed data. The function computes the following indices at each draw
\\i\\:

- **BRMSEA** (Bayesian Root Mean Square Error of Approximation):
  \$\$\text{BRMSEA}\_i = \sqrt{\max \[0, (D_i^{obs} - p^\*) / (df_B
  \cdot N)\]}\$\$

- **BGamma** (\\\hat{\Gamma}\\) and **Adjusted BGamma**
  (\\\hat{\Gamma}\_{adj}\\): Bayesian analogs of the GFI.
  \$\$\hat{\Gamma}\_i = J / \[J + (2 / N) \cdot (D_i^{obs} - p^\*)\]\$\$
  \$\$\hat{\Gamma}\_{adj,i} = 1 - (p^\* / df_B) \cdot (1 -
  \hat{\Gamma}\_i)\$\$

- **BMc** (Bayesian McDonald's Centrality Index): \$\$\text{BMc}\_i =
  \exp \[-(1 / (2 \cdot N)) \cdot (D_i^{obs} - p^\*)\]\$\$

- **Chi2_ppp** (Posterior Predictive P-value): The proportion of MCMC
  iterations where the deviance of replicated data exceeds the observed
  deviance \\P(D^{rep} \geq D^{obs})\\. Values closer to 0.50 indicate
  excellent fit, while values \\\< .05\\ suggest poor fit.

- **SRMR**: The Standardized Root Mean Square Residual, computed
  explicitly at each MCMC draw.

**Incremental Fit Indices**

Incremental indices compare the hypothesized model (\\M_H\\) against a
worst-case baseline independence model (\\M_0\\), where all observed
variables are assumed to be uncorrelated. *Note: The function
automatically estimates this null model under the hood using the same
Stan settings, which may increase computation time.*

- **BCFI** (Bayesian Comparative Fit Index): \$\$\text{BCFI}\_i = 1 -
  \[(D\_{H,i}^{obs} - p^\*) / (D\_{0,i}^{obs} - p^\*)\]\$\$

- **BTLI** (Bayesian Tucker-Lewis Index / Non-Normed Fit Index):
  \$\$\text{BTLI}\_i = \[ (D\_{0,i}^{obs} - p^\*) / df\_{B0} -
  (D\_{H,i}^{obs} - p^\*) / df\_{BH} \] / \[ (D\_{0,i}^{obs} - p^\*) /
  df\_{B0} - 1 \]\$\$

*(Note: Both BCFI and BTLI are strictly bounded to the \\\[0, 1\]\\
interval in the final output).*

**Predictive Fit Indices**

Out-of-sample predictive accuracy is assessed using PSIS-LOO
cross-validation:

- **ELPD**: Expected Log Pointwise Predictive Density.

- **LOOIC**: Leave-One-Out Information Criterion (\\-2 \cdot
  \text{ELPD}\\).

- **p_loo**: Estimated effective number of parameters (\\pD\\).

## References

Garnier-Villarreal, M., & Jorgensen, T. D. (2020). Adapting Fit Indices
for Bayesian Structural Equation Modeling: Comparison to Maximum
Likelihood. *Psychological Methods, 25*(1), 46-70.
<https://doi.org/10.1037/met0000224>

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing, 27*(5), 1413-1432.
<https://doi.org/10.1007/s11222-016-9696-4>

## Examples

``` r
if (FALSE) { # \dontrun{
# --------------------------------------------------------------- #
#    1. Bayesian Fit Measures after fitting Bayesian EFA model    #
# --------------------------------------------------------------- #

# Fit Bayesian EFA model to the famous Grant-White School Data (Holzinger y Swineford , 1939)
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

# Compute Bayesian Fit Measures from a fitted befa object
bayesian_fit_measures <- befa_fit_measures(befa_fit)

# Posterios summaries (NAs when posterior draws have sds close to zero)
bayesian_fit_measures$summary

# # A tibble: 9 × 10
#   variable       mean   median      sd     mad       q5      q95  rhat ess_bulk ess_tail
#   <chr>         <dbl>    <dbl>   <dbl>   <dbl>    <dbl>    <dbl> <dbl>    <dbl>    <dbl>
# 1 Chi2        47.8     47.1    7.24    7.19     37.2     60.9     1.00    1276.    2504.
# 2 Chi2_Null  919.     919.     0       0       919.     919.     NA         NA       NA
# 3 BRMSEA       0.0576   0.0594 0.0209  0.0195    0.0196   0.0890  1.00    1303.    2504.
# 4 BGamma       0.991    0.992  0.00517 0.00522   0.982    0.999   1.00    1303.    2504.
# 5 Adj_BGamma   0.970    0.972  0.0179  0.0180    0.938    0.997   1.00    1303.    2504.
# 6 BMc          0.981    0.982  0.0116  0.0117    0.960    0.998   1.00    1303.    2504.
# 7 SRMR         0.0465   0.0456 0.00809 0.00773   0.0349   0.0612  1.00    1836.    2978.
# 8 BCFI         0.987    0.987  0.00810 0.00814   0.972    0.999   1.00    1303.    2504.
# 9 BTLI         0.986    0.997  0.0208  0.00381   0.941    1       1.00    1466.      NA

# Posterior draws of fit measures (transform to matrix)
posterior::as_draws_matrix(bayesian_fit_measures$fit_draws)

# # A draws_matrix: 1000 iterations, 4 chains, and 9 variables
#     variable
# draw Chi2 Chi2_Null BRMSEA BGamma Adj_BGamma  BMc  SRMR BCFI
#   1    39       919  0.031   1.00       0.99 1.00 0.038 1.00
#   2    52       919  0.070   0.99       0.96 0.97 0.049 0.98
#   3    43       919  0.048   0.99       0.98 0.99 0.041 0.99
#   4    60       919  0.088   0.98       0.94 0.96 0.051 0.97
#   5    50       919  0.067   0.99       0.96 0.98 0.049 0.98
#   6    49       919  0.063   0.99       0.97 0.98 0.044 0.99
#   7    53       919  0.074   0.99       0.96 0.97 0.048 0.98
#   8    57       919  0.082   0.98       0.95 0.97 0.047 0.98
#   9    60       919  0.088   0.98       0.94 0.96 0.056 0.97
#   10   61       919  0.090   0.98       0.94 0.96 0.049 0.97
# # ... with 3990 more draws, and 1 more variables

# We can also plot the posterior distribution of fit measures
hist(posterior::as_draws_matrix(bayesian_fit_measures$fit_draws)[, "SRMR"])

# PSIS-LOO Expected log-predictive density value, uncertainty, and pareto-k values
bayesian_fit_measures$loo_object

# Computed from 4000 by 301 log-likelihood matrix.
#
#          Estimate   SE
# elpd_loo  -3416.9 42.5
# p_loo        25.6  1.8
# looic      6833.9 85.0
# ------
# MCSE of elpd_loo is 0.1.
# MCSE and ESS estimates assume independent draws (r_eff=1).
#
# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

# -----------------------------------------------------------------------------

# -------------------------------------------------------------- #
#    2. Bayesian Fit Measures when fitting Bayesian EFA model    #
# -------------------------------------------------------------- #

# Fit Bayesian EFA model to the famous Grant-White School Data (Holzinger y Swineford , 1939)
befa_fit <- befa(
  data = HS_data,
  n_factors = 3,
  factor_scores = FALSE,
  compute_fit_indices = TRUE,
  compute_reliability = FALSE,
  backend = "rstan",
  seed = 17,
  chains = 4,
  parallel_chains = 4
)

# befa_fit_measures output is inside the "fit_indices" object
# befa_fit$fit_indices

# ----------------------------------------------------------------------------
} # }
```
