# Efficient Rotation-Sign-Permutation (E-RSP) Alignment

Applies the Efficient Rotation-Sign-Permutation (E-RSP) algorithm to
MCMC draws of factor loadings. This post-processing step resolves
rotational indeterminacy, label switching, and sign reflections across
posterior draws, enabling the computation of valid and interpretable
posterior summaries.

## Usage

``` r
rsp_align(
  lambda_draws,
  n_items,
  n_factors,
  n_chains = 1,
  format = NULL,
  max_iter = 1000,
  threshold = 1e-06,
  add_names = TRUE
)
```

## Arguments

- lambda_draws:

  A numeric 2D matrix of dimensions S \\\times\\ (J \\\times\\ M), where
  S is the number of posterior draws. Each row must contain a flattened
  loading matrix from a single MCMC iteration.

- n_items:

  Integer. The number of observed variables (J).

- n_factors:

  Integer. The number of latent factors (M).

- n_chains:

  Integer. The number of MCMC chains used for the posterior draws. This
  argument is required to correctly handle the structure of the draws.
  The function assumes that the `lambda_draws` matrix is constructed by
  stacking the draws of each chain row-by-row. That is, the first
  \\S/n\_{chains}\\ rows correspond to the first chain, the next
  \\S/n\_{chains}\\ to the second chain, and so on.

- format:

  Character string specifying how the loading matrix was flattened into
  a row vector for each MCMC draw (i.e., the order of the values in each
  row of `lambda_draws`). The aligned output will be returned in this
  same format. **This argument is required** (no default is provided).
  See the **Flattening Formats** section in Details for a visual
  explanation. Options are:

  - `"column_major"`: Elements are ordered column by column (factor by
    factor), corresponding to `as.vector(Lambda)` in R.

  - `"row_major"`: Elements are ordered row by row (item by item),
    corresponding to `as.vector(t(Lambda))` in R.

- max_iter:

  Integer. The maximum number of iterations for the RSP algorithm to
  reach convergence.

- threshold:

  Numeric. The convergence threshold for the alignment objective
  function (Frobenius discrepancy).

- add_names:

  Logical. If `TRUE` (default), assigns Stan-format column names (e.g.,
  `Lambda[1,1]`, `Lambda[2,1]`) to the output aligned matrix. If
  `FALSE`, the function checks whether `lambda_draws` already has column
  names: if it does, those names are preserved in the output; if it does
  not, `add_names` is automatically set to `TRUE`.

## Value

A list containing:

- `Lambda_hat_mcmc`: A numeric matrix of the aligned loading draws,
  maintaining the exact dimensions and format of the input
  `lambda_draws`. The attribute `nchains` is added to this matrix,
  indicating the number of chains (as passed in `n_chains`).

- `Lambda_star`: A \\J \times M\\ numeric matrix representing the final
  reference configuration (the posterior mean of the aligned draws).

- `objective`: A numeric value indicating the final alignment objective
  (the total Frobenius discrepancy across all draws).

- `sign_vectors`: A matrix of dimensions \\S \times M\\ containing the
  sign flips applied to each draw.

- `perm_vectors`: A matrix of dimensions \\S \times M\\ containing the
  permutations applied to each draw.

- `summary`: A `draws_summary` data frame computed via
  [`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html).

## Details

**Why Alignment is Necessary**

Factor models are invariant to orthogonal transformations: if
\\\mathbf{\Lambda}\\ is a valid loading matrix, then
\\\mathbf{\Lambda}^\bullet = \mathbf{\Lambda} \cdot \mathbf{Q}\\ yields
an identical model-implied covariance matrix for any orthogonal matrix
\\\mathbf{Q}\\ where \\\mathbf{Q} \cdot \mathbf{Q}^\top =
\mathbf{I}\_M\\. During MCMC sampling, this invariance creates a
multimodal posterior where the sampler explores equivalent, yet
incompatible, rotational orientations. Consequently, naive posterior
summaries (such as averaging across draws) mix these modes, yielding
uninterpretable loading estimates that artificially cancel out toward
zero.

**The Efficient RSP Algorithm (E-RSP)**

This function implements the Efficient Rotation-Sign-Permutation (E-RSP)
algorithm (Rey-Sáez et al., 2026), an optimized version of the exact RSP
method originally proposed by Papastamoulis and Ntzoufras (2022). While
the original approach becomes computationally prohibitive as the number
of latent factors increases, E-RSP reduces the alignment task to a
Linear Assignment Problem (LAP). This guarantees a globally optimal
solution that scales to high-dimensional models with negligible
computational cost.

The algorithm proceeds in two stages:

1.  **Continuous Alignment (Varimax)**: Each posterior draw is rotated
    to a canonical simple structure to fix the continuous rotational
    degree of freedom.

2.  **Discrete Alignment (Signed-Permutation)**: The remaining discrete
    ambiguity (column permutations and sign flips) is resolved by
    mapping each draw to a common reference configuration.

**Flattening Formats**

The `lambda_draws` input must be a 2D matrix of dimensions S \\\times\\
(J \\\times\\ M), where each row represents a flattened loading matrix
from a single MCMC draw. The `format` argument simply tells the function
how your MCMC software flattened the original \\J \times M\\ loading
matrix \\\mathbf{\Lambda}\\ into that row vector. The aligned output
will be returned in this exact same flattened format.

To illustrate, consider a \\4 \times 3\\ loading matrix (4 items, 3
factors): \$\$\mathbf{\Lambda} = \begin{bmatrix} \lambda\_{11} &
\lambda\_{12} & \lambda\_{13} \\ \lambda\_{21} & \lambda\_{22} &
\lambda\_{23} \\ \lambda\_{31} & \lambda\_{32} & \lambda\_{33} \\
\lambda\_{41} & \lambda\_{42} & \lambda\_{43} \end{bmatrix}\$\$

Depending on the `format` specified, the input rows must be structured
as follows:

- **Column-major** (`format = "column_major"`): Elements are filled
  column by column (factor by factor). This corresponds to the
  `as.vector(Lambda)` behavior in R and is the default output layout for
  Stan: \$\$\text{vec}(\mathbf{\Lambda}) = \big\[
  \underbrace{\lambda\_{11}, \lambda\_{21}, \lambda\_{31},
  \lambda\_{41}}\_{\text{Factor 1}}, \\ \underbrace{\lambda\_{12},
  \lambda\_{22}, \lambda\_{32}, \lambda\_{42}}\_{\text{Factor 2}}, \\
  \underbrace{\lambda\_{13}, \lambda\_{23}, \lambda\_{33},
  \lambda\_{43}}\_{\text{Factor 3}} \big\]\$\$

- **Row-major** (`format = "row_major"`): Elements are filled row by row
  (item by item). This corresponds to `as.vector(t(Lambda))` in R:
  \$\$\text{vec}(\mathbf{\Lambda}^\top) = \big\[
  \underbrace{\lambda\_{11}, \lambda\_{12}, \lambda\_{13}}\_{\text{Item
  1}}, \\ \underbrace{\lambda\_{21}, \lambda\_{22},
  \lambda\_{23}}\_{\text{Item 2}}, \\ \underbrace{\lambda\_{31},
  \lambda\_{32}, \lambda\_{33}}\_{\text{Item 3}}, \\
  \underbrace{\lambda\_{41}, \lambda\_{42}, \lambda\_{43}}\_{\text{Item
  4}} \big\]\$\$

## References

Papastamoulis, P., & Ntzoufras, I. (2022). On the identifiability of
Bayesian factor analytic models. *Statistics and Computing, 32*(2), 23.
<https://doi.org/10.1007/s11222-022-10084-4>

Rey-Sáez, R. & Revuelta, J. (2026). *An Efficient
Rotation-Sign-Permutation Algorithm to Solve Rotational Indeterminacy in
Bayesian Exploratory Factor Analysis*. Manuscript submitted for
publication. <https://osf.io/5dutv/>

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit Bayesian EFA model
befa_fit <- befa(
  data = HS_data,
  n_factors = 3,
  rotate = "none",
  factor_scores = FALSE,
  compute_fit_indices = FALSE,
  compute_reliability = FALSE,
  backend = "rstan",
  seed = 17,
  chains = 4,
  parallel_chains = 4
)

# Extract unrotated posterior draws
lambda_unrotated <- extract_posterior_draws(befa_fit, pars = "Lambda")

# See multimodality due to rotational indeterminacy
hist(lambda_unrotated[, 1],
  breaks = 100, col = "steelblue2",
  main = "Rotation indeterminacy", xlab = "Lambda[1,1]"
)

# See that columns are ordered following a column-major order
# see Details on rsp_align function
colnames(lambda_unrotated)

#  [1] "Lambda[1,1]" "Lambda[2,1]" "Lambda[3,1]" "Lambda[4,1]" "Lambda[5,1]" "Lambda[6,1]" "Lambda[7,1]"
#  [8] "Lambda[8,1]" "Lambda[9,1]" "Lambda[1,2]" "Lambda[2,2]" "Lambda[3,2]" "Lambda[4,2]" "Lambda[5,2]"
# [15] "Lambda[6,2]" "Lambda[7,2]" "Lambda[8,2]" "Lambda[9,2]" "Lambda[1,3]" "Lambda[2,3]" "Lambda[3,3]"
# [22] "Lambda[4,3]" "Lambda[5,3]" "Lambda[6,3]" "Lambda[7,3]" "Lambda[8,3]" "Lambda[9,3]"

# Now, let's align posterior draws
lambda_aligned <- rsp_align(
  lambda_draws = lambda_unrotated,
  n_items = ncol(HS_data),
  n_factors = 3,
  n_chains = 4,
  format = "column_major"
)

# Now, rotational indeterminacy has been solved!
hist(lambda_aligned$Lambda_hat_mcmc[, 1],
  breaks = 100, col = "steelblue2",
  main = "Aligned posterior distribution", xlab = "Lambda[1,1]"
)
} # }
```
