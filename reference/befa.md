# Bayesian Exploratory Factor Analysis

Estimates an Exploratory Factor Analysis model using Hamiltonian Monte
Carlo (Stan). Posterior samples are post-processed via Varimax rotation
and alignment (RSP algorithm) to resolve the rotational indeterminacy
inherent to factor models.

## Usage

``` r
befa(
  data = NULL,
  n_factors = NULL,
  ordered = FALSE,
  sample_nobs = NULL,
  sample_mean = NULL,
  sample_cov = NULL,
  sample_cor = NULL,
  model = c("cor", "cov", "raw"),
  lambda_prior = c("unit_vector", "normal"),
  missing = c("listwise", "FIML"),
  rotate = c("varimax", "none"),
  rsp_args = NULL,
  loo_args = NULL,
  factor_scores = FALSE,
  backend = c("rstan", "cmdstanr"),
  prior = list(),
  compute_fit_indices = TRUE,
  compute_reliability = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  Numeric matrix or data.frame (N rows, J columns). Raw observations.
  Mutually exclusive with `sample_cov`/`sample_cor`.

- n_factors:

  Integer. Number of latent factors to extract.

- ordered:

  Logical. If TRUE, fits an ordinal model. Currently not implemented.

- sample_nobs:

  Integer. Sample size. Required when using summary statistics instead
  of `data`.

- sample_mean:

  Numeric vector (length J). Sample means. Used with `model="raw"` and
  summary stats.

- sample_cov:

  Numeric matrix (J x J). Sample covariance matrix. Mutually exclusive
  with `data`.

- sample_cor:

  Numeric matrix (J x J). Sample correlation matrix. Mutually exclusive
  with `data`.

- model:

  Character. Model parameterization:

  - `"cor"` (default): Estimates the model-implied correlation matrix.

  - `"cov"`: Estimates the model-implied covariance matrix.

  - `"raw"`: Full model. Estimates the model-implied mean-vector and
    covariance matrix.

- lambda_prior:

  Character. Prior family for loadings:

  - `"unit_vector"` (default): Unit-vector prior distribution proposed
    by Rey-Sáez et al., 2025 (see details)

  - `"normal"`: Independent normal priors for each factor loading.
    Normal priors are not available when `model="cor"`.

- missing:

  Character. Method for handling missing data:

  - `"listwise"` (default): Complete-case analysis (removes rows with
    any NA).

  - `"FIML"`: Full Information Maximum Likelihood. Uses all available
    data, estimating parameters using the observed data pattern for each
    case. Requires raw data.

- rotate:

  Character. Rotation criterion for post-processing:

  - `"varimax"` (default): Apply Varimax rotation and RSP alignment.

  - `"none"`: No rotation or alignment applied.

- rsp_args:

  List. Arguments for the efficient RSP alignment algorithm:

  - `max_iter`: Maximum iterations (default 1000).

  - `threshold`: Convergence threshold (default 1e-6).

- loo_args:

  List. Arguments passed to
  [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) for
  leave-one-out cross-validation when `compute_fit_indices = TRUE`.
  Common arguments include:

  - `cores`: Number of cores to use for parallelization.

  - `r_eff`: Logical. If TRUE, computes relative effective sample sizes.

- factor_scores:

  Logical. If TRUE, computes and stores factor scores. Requires `data`.

- backend:

  Character. Stan interface: `"rstan"` (default) or `"cmdstanr"`.

- prior:

  Named list. User-defined prior hyperparameters (override defaults):

  - `xi`: Repulsion parameter (unit_vector prior). Default: `100`. **Not
    recommended to be modified by the user**.

  - `nu`: c(mean, sd) for intercepts (raw model). Default: `c(0, 40)`
    -\> N(0, 40).

  - `sigma`: c(df, loc, scale) for residual SDs (cov/raw models).
    Default: `c(3, 0, 2.5)` -\> half-t(3, 0, 2.5).

  - `h2`: c(alpha, beta) for communalities (unit_vector prior). Default:
    `c(1, 1)` -\> Beta(1, 1). Using this default setting implies a
    **uniform joint distribution over the factor loadings** when
    `lambda_prior = "unit_vector"`.

  - `lambda`: c(mean, sd) for loadings (normal prior). Default:
    `c(0, 10)` -\> N(0, 10). For the unidimensional model, c(alpha,
    beta) -\> Beta(1, 1).

  - `psi`: c(alpha, beta) for uniquenesses (normal prior). Default:
    `c(0.5, 0.5)` -\> InvGamma(0.5, 0.5).

- compute_fit_indices:

  Logical. If TRUE (default), computes Bayesian fit indices after
  estimation.

- compute_reliability:

  Logical. If TRUE (default), computes Omega reliability coefficients.

- verbose:

  Logical. If TRUE (default), prints progress messages.

- ...:

  Additional arguments passed to Stan sampler (iter, chains, warmup,
  seed, etc.).

## Value

An object of class `befa` containing:

- `stanfit`: The Stan fit object with aligned posterior samples.

- `stan_data`: Data list passed to Stan.

- `n_factors`, `model_type`, `lambda_prior`, `rotation`: Model
  specifications.

- `missing`: Missing data method used (`"listwise"` or `"FIML"`).

- `has_missing`: Logical. Whether the data contained missing values.

- `priors_used`: Resolved prior hyperparameters.

- `options`: List with `factor_scores`, `ordered`, `rsp_args`, and
  `loo_args`.

- `rsp_objective`: Final RSP alignment objective value (NA if
  `rotate = "none"`).

- `call`: The matched function call.

- `reliability`: Omega coefficients (if `compute_reliability = TRUE`).

- `fit_indices`: Bayesian fit measures (if
  `compute_fit_indices = TRUE`).

## Details

***The Factor Model***

The function estimates a linear factor model where each observed vector
\\\mathbf{y}\_i\\ is decomposed as: \$\$\mathbf{y}\_i =
\boldsymbol{\nu} + \mathbf{\Lambda} \cdot \boldsymbol{\eta}\_i +
\boldsymbol{\epsilon}\_i\$\$

where:

- \\\mathbf{y}\_i\\ is a \\J\\-dimensional observed response vector.

- \\\boldsymbol{\nu}\\ is a \\J \times 1\\ vector of intercepts.

- \\\mathbf{\Lambda}\\ is the \\J \times M\\ matrix of factor loadings.

- \\\boldsymbol{\eta}\_i \sim \mathcal{N}(\mathbf{0}, \mathbf{I}\_M)\\
  represents the latent factors.

- \\\boldsymbol{\epsilon}\_i \sim \mathcal{N}(\mathbf{0},
  \mathbf{\Psi})\\ represents the idiosyncratic residuals
  (uniquenesses), with \\\mathbf{\Psi}\\ being a diagonal matrix.

Under the assumption of no missing data (or using listwise deletion),
the marginal distribution of the observed data is: \$\$\mathbf{y}\_i
\mid \boldsymbol{\nu}, \mathbf{\Lambda}, \mathbf{\Psi} \sim
\mathcal{MVN}(\boldsymbol{\nu}, \mathbf{\Sigma})\$\$

where \\\mathbf{\Sigma} = \mathbf{\Lambda}\cdot\mathbf{\Lambda}^\top +
\mathbf{\Psi}\\ is the model-implied covariance matrix.

***Model Parameterizations***

The `befa` function allows users to fit different versions of the
marginal model through the following parameterizations:

- **Correlation Model** (`model = "cor"`): This parameterization assumes
  the observed data are standardized; thus, the intercept vector
  \\\boldsymbol{\nu}\\ is not estimated. The model-implied covariance
  matrix \\\mathbf{\Sigma}\\ simplifies to the model-implied correlation
  matrix \\\mathbf{R}\\, decomposed as: \$\$\mathbf{R} =
  \mathbf{\Lambda} \cdot \mathbf{\Lambda}^\top + \mathbf{\Psi}\$\$ where
  \\\mathbf{\Lambda}\\ contains the standardized factor loadings and
  \\\mathbf{\Psi}\\ represents the proportion of unique variance
  (uniquenesses). Since the elements of \\\mathbf{\Lambda}\\ are
  constrained to the range \\\[-1, 1\]\\ and are mutually dependent (as
  the diagonal elements of \\\mathbf{R}\\ must sum to 1),
  `lambda_prior = "normal"` is not available. Instead, the unit-vector
  prior (`lambda_prior = "unit_vector"`) is employed to ensure a valid
  joint distribution by accounting for this dependency.

- **Covariance Model** (`model = "cov"`): This model focuses on the
  covariance structure while assuming the intercept vector
  \\\boldsymbol{\nu}\\ is zero. It can be implemented in two ways:

  - *Standardized Structure with Scale*: When using
    `lambda_prior = "unit_vector"`, the model introduces a vector of
    marginal standard deviations \\\boldsymbol{\sigma}\\. The covariance
    matrix is reconstructed by scaling the standardized structure:
    \$\$\mathbf{\Sigma} = \mathbf{D}\_{\sigma} \cdot (\mathbf{\Lambda}
    \cdot \mathbf{\Lambda}^\top + \mathbf{\Psi}) \cdot
    \mathbf{D}\_{\sigma}\$\$ where \\\mathbf{D}\_{\sigma} =
    \text{diag}(\boldsymbol{\sigma})\\ and \\\mathbf{\Lambda}\\ remains
    in the correlation metric.

  - *Unstandardized (Metric) Structure*: When using
    `lambda_prior = "normal"`, the loadings are estimated directly in
    the scale of the observed variables. The scale is embedded within
    \\\mathbf{\Lambda}\\, and the covariance matrix follows the
    classical fundamental equation: \$\$\mathbf{\Sigma} =
    \mathbf{\Lambda} \cdot \mathbf{\Lambda}^\top + \mathbf{\Psi}\$\$
    where \\\mathbf{\Psi}\\ represents uniquenesses in the original
    variance scale.

- **Full Unstandardized Model** (`model = "raw"`): Unlike the previous
  parameterizations, this model estimates the full marginal distribution
  of the observed data, which includes estimating the intercept vector
  \\\boldsymbol{\nu}\\ (i.e., the item means). The covariance structure
  \\\mathbf{\Sigma}\\ is modeled exactly as in the Covariance Model
  (`model = "cov"`). Depending on the `lambda_prior` selected, the
  covariance matrix is reconstructed either by scaling a standardized
  factor structure (`lambda_prior = "unit_vector"`) or by estimating the
  loadings directly in the original metric (`lambda_prior = "normal"`).

***Prior distributions***

The following arguments allow users to specify the prior distributions
for the model parameters. The exact parameters estimated depend on the
chosen `model` and `lambda_prior`.

- **Unit-Vector Prior for Loadings** (`lambda_prior = "unit_vector"`):
  Rather than placing priors directly on the individual loadings, the
  user only needs to specify a prior on the communality, \\h_j^2 \sim
  \text{Beta}(\alpha, \beta)\\, which represents the proportion of
  variance explained by the factors. This prior is coherently translated
  to the factor loadings via the decomposition \\\boldsymbol{\lambda}\_j
  = \sqrt{h_j^2} \cdot \mathbf{z}\_j\\, where \\\mathbf{z}\_j\\ is a
  vector uniformly distributed on the unit hypersphere. **By default
  (using a \\\text{Beta}(1, 1)\\ prior on the communalities), this
  approach implicitly assumes a uniform joint distribution over the
  valid space of the factor loadings.** This accounts for the dependency
  among loadings and strictly prevents Heywood cases (negative residual
  variances). See Rey-Sáez et al. (2025) for details.

- **Normal Prior for Loadings** (`lambda_prior = "normal"`): This
  approach assigns standard independent normal priors directly to each
  loading: \\\lambda\_{jm} \sim \mathcal{N}(\mu\_\lambda,
  \sigma\_\lambda)\\. While intuitive and commonly used, it does not
  constrain the variance partitioning and is unavailable for the
  Correlation Model (`model = "cor"`).

- **Prior for Uniquenesses** (\\\mathbf{\Psi}\\): This prior is only
  required when using `lambda_prior = "normal"` in the Covariance
  (`"cov"`) or Full Unstandardized (`"raw"`) models. It defines the
  prior distribution for the residual variances, \\\psi\_{jj}\\, which
  follow an inverse-gamma distribution: \\\psi\_{jj} \sim
  \mathcal{IG}(\alpha, \beta)\\. *(Note: When using
  `lambda_prior = "unit_vector"`, this prior is ignored because the
  uniquenesses in the standardized structure are deterministically
  defined by the communality as \\1 - h_j^2\\).*

- **Prior for Scale Parameters** (\\\boldsymbol{\sigma}\\): This prior
  is only used when scaling the standardized structure back to the
  original metric (i.e., using `lambda_prior = "unit_vector"` with
  `model = "cov"` or `model = "raw"`). It assigns a Student-t prior
  distribution to the marginal standard deviations of each observed
  variable, \\\sigma_j\\: \\\sigma_j \sim \text{Student-t}(\nu\_\sigma,
  \mu\_\sigma, \sigma\_\sigma)\\.

- **Prior for Intercepts** (\\\boldsymbol{\nu}\\): This prior is
  exclusively used in the Full Unstandardized Model (`model = "raw"`),
  regardless of the chosen `lambda_prior`. It assigns a normal prior
  distribution to the item means, \\\nu_j\\: \\\nu_j \sim
  \mathcal{N}(\mu\_\nu, \sigma\_\nu)\\.

## References

Garnier-Villarreal, M., & Jorgensen, T. D. (2020). Adapting fit indices
for Bayesian structural equation modeling: Comparison to maximum
likelihood. *Psychological methods, 25*(1), 46–70.
<https://doi.org/10.1037/met0000224>

Holzinger, K. J., and F. A. Swineford. 1939. *A Study of Factor
Analysis: The Stability of a Bi-Factor Solution.* Supplementary
Educational Monograph 48. Chicago: University of Chicago Press.

Rey-Sáez, R., Franco-Martínez, A., Revuelta, J., & Vadillo, M. A.
(2025). *A Unified Framework for Psychometrics in Experimental
Psychology: The Standardized Generalized Hierarchical Factor Model*.
PsyArXiv. <https://doi.org/10.31234/osf.io/gv6k7_v1>

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit Bayesian EFA model to the famous Grant-White School Data (Holzinger & Swineford, 1939)
befa_fit <- befa(
  data = HS_data,
  n_factors = 3,
  model = "cor",
  lambda_prior = "unit_vector",
  rotate = "varimax",
  factor_scores = TRUE,
  compute_fit_indices = TRUE,
  compute_reliability = TRUE,
  backend = "rstan",
  iter_sampling = 1000,
  iter_warmup = 1000,
  chains = 4,
  parallel_chains = 4,
  seed = 17
)

# Posterior draws of model parameters
extract_posterior_draws(befa_fit, pars = "Lambda", format = "matrix")

# # A draws_matrix: 1000 iterations, 4 chains, and 27 variables
#     variable
# draw Lambda[1,1] Lambda[2,1] Lambda[3,1] Lambda[4,1] Lambda[5,1] Lambda[6,1]
#   1         0.63        0.50        0.66       0.053     0.02543        0.19
#   2         0.62        0.48        0.68       0.132     0.01635        0.11
#   3         0.55        0.43        0.69       0.199     0.00023        0.17
#   4         0.44        0.42        0.55       0.096    -0.00680        0.14
#   5         0.59        0.51        0.69       0.143    -0.00072        0.22
#   6         0.65        0.42        0.57       0.075     0.03311        0.17
#   7         0.58        0.49        0.74       0.106    -0.00305        0.20
#   8         0.69        0.42        0.73       0.084     0.07283        0.19
#   9         0.45        0.45        0.64       0.115     0.00454        0.15
#   10        0.61        0.49        0.64       0.093     0.06957        0.19
# # ... with 3990 more draws, and 21 more variables

# Fast summaries using rvars format
extract_posterior_draws(befa_fit, pars = "Lambda", format = "rvars")

# # A draws_rvars: 1000 iterations, 4 chains, and 1 variables
# $Lambda: rvar<1000,4>[9,3] mean ± sd:
#       [,1]            [,2]            [,3]
#  [1,]  0.594 ± 0.063   0.313 ± 0.048   0.122 ± 0.060
#  [2,]  0.459 ± 0.062   0.130 ± 0.054  -0.038 ± 0.067
#  [3,]  0.648 ± 0.063   0.076 ± 0.048   0.108 ± 0.054
#  [4,]  0.110 ± 0.039   0.833 ± 0.024   0.073 ± 0.040
#  [5,]  0.029 ± 0.037   0.862 ± 0.024   0.068 ± 0.038
#  [6,]  0.156 ± 0.040   0.809 ± 0.025   0.062 ± 0.041
#  [7,] -0.053 ± 0.050   0.099 ± 0.048   0.681 ± 0.077
#  [8,]  0.171 ± 0.065   0.075 ± 0.046   0.693 ± 0.074
#  [9,]  0.401 ± 0.071   0.165 ± 0.049   0.492 ± 0.067

# Full posterior summaries
posterior_summaries(befa_fit, pars = "h2")

# # A tibble: 9 × 10
#   variable  mean median     sd    mad    q5   q95  rhat ess_bulk ess_tail
#   <chr>    <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
# 1 h2[1]    0.476  0.474 0.0645 0.0622 0.370 0.584 1.00     4808.    2810.
# 2 h2[2]    0.240  0.237 0.0569 0.0568 0.148 0.338 1.00     4910.    2645.
# 3 h2[3]    0.447  0.442 0.0791 0.0770 0.325 0.582 1.00     4352.    2399.
# 4 h2[4]    0.715  0.717 0.0364 0.0356 0.654 0.774 1.00     5036.    2427.
# 5 h2[5]    0.751  0.752 0.0386 0.0384 0.686 0.812 1.00     5169.    3079.
# 6 h2[6]    0.687  0.688 0.0358 0.0360 0.626 0.743 1.00     4945.    2674.
# 7 h2[7]    0.486  0.479 0.103  0.0990 0.332 0.669 1.00     3504.    2121.
# 8 h2[8]    0.527  0.523 0.0873 0.0846 0.390 0.674 1.00     3800.    2562.
# 9 h2[9]    0.442  0.441 0.0560 0.0549 0.350 0.533 1.000    4101.    2681.

# # Model summaries (inspired in psych package output)
# summary(befa_fit, cutoff = 0.3, signif_stars = TRUE)
#
# Table 1. Factor Loadings (Pattern Matrix)
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Variable       F1     F2     F3    h2    u2  Rhat  EssBulk  EssTail
# ———————————————————————————————————————————————————————————————————
# Item_1      0.59*  0.31*         0.47  0.53  1.00     3458     2849
# Item_2      0.46*                0.23  0.77  1.00     4053     2741
# Item_3      0.65*                0.44  0.56  1.00     4146     2274
# Item_4             0.83*         0.71  0.29  1.00     4156     2607
# Item_5             0.86*         0.75  0.25  1.00     4179     2983
# Item_6             0.81*         0.68  0.32  1.00     4291     2817
# Item_7                    0.68*  0.48  0.52  1.00     3343     2081
# Item_8                    0.69*  0.52  0.48  1.00     3614     2593
# Item_9      0.40*         0.49*  0.43  0.57  1.00     3014     2757
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Note: varimax rotation applied. Diagnostics show worst-case values
# across factors (max Rhat, min ESS). The 3 latent factors accounted
# for 52.2% of total variance. (*) 95% Credible Interval excludes 0.
# Loadings with absolute values < 0.30 are hidden.
#
# Table 2. Bayesian Fit Measures
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Index         Estimate     SD    CI_Low   CI_High
# —————————————————————————————————————————————————
# Chi2             47.75   7.24     35.81     63.71
# Chi2_ppp          0.12
# Chi2_Null       918.85   0.00    918.85    918.85
# BRMSEA            0.06   0.02      0.00      0.09
# BGamma            0.99   0.01      0.98      1.00
# Adj_BGamma        0.97   0.02      0.93      1.00
# BMc               0.98   0.01      0.96      1.00
# SRMR              0.05   0.01      0.03      0.06
# BCFI              0.99   0.01      0.97      1.00
# BTLI              0.99   0.02      0.93      1.00
# ELPD          -3416.93  42.49  -3500.21  -3333.64
# LOOIC          6833.85  84.98   6667.29   7000.42
# p_loo            25.58   1.82     22.01     29.14
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Note: Intervals are 95% Credible Intervals. PPP:
# Posterior Predictive p-value (Ideal > .05).
# p_loo/LOOIC derived from PSIS-LOO.
#
# Table 3. Factor Reliability (Coefficient Omega)
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Factor    Estimate    SD  CI_Low  CI_High
# —————————————————————————————————————————
# F1            0.60  0.04    0.52     0.67
# F2            0.73  0.02    0.69     0.76
# F3            0.55  0.04    0.45     0.62
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Full Scale Omega Total: 0.84 [0.82,
# 0.86]. Omega coefficients use the full
# posterior distribution.
} # }
```
