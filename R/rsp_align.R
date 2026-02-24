#' Efficient Rotation-Sign-Permutation (E-RSP) Alignment
#'
#' Applies the Efficient Rotation-Sign-Permutation (E-RSP) algorithm to MCMC draws
#' of factor loadings. This post-processing step resolves rotational indeterminacy,
#' label switching, and sign reflections across posterior draws, enabling the
#' computation of valid and interpretable posterior summaries.
#'
#' @param lambda_draws A numeric 2D matrix of dimensions S \eqn{\times} (J \eqn{\times} M),
#'   where S is the number of posterior draws. Each row must contain a flattened
#'   loading matrix from a single MCMC iteration.
#' @param n_items Integer. The number of observed variables (J).
#' @param n_factors Integer. The number of latent factors (M).
#' @param n_chains Integer. The number of MCMC chains used for the posterior draws.
#'   This argument is required to correctly handle the structure of the draws.
#'   The function assumes that the `lambda_draws` matrix is constructed by stacking
#'   the draws of each chain row-by-row. That is, the first \eqn{S/n_{chains}}
#'   rows correspond to the first chain, the next \eqn{S/n_{chains}} to the
#'   second chain, and so on.
#' @param format Character string specifying how the loading matrix was flattened
#'   into a row vector for each MCMC draw (i.e., the order of the values in each
#'   row of `lambda_draws`). The aligned output will be returned in this same
#'   format. **This argument is required** (no default is provided). See the
#'   **Flattening Formats** section in Details for a visual explanation. Options are:
#'   * `"column_major"`: Elements are ordered column by column (factor by factor),
#'     corresponding to `as.vector(Lambda)` in R.
#'   * `"row_major"`: Elements are ordered row by row (item by item),
#'     corresponding to `as.vector(t(Lambda))` in R.
#' @param max_iter Integer. The maximum number of iterations for the RSP algorithm
#'   to reach convergence.
#' @param threshold Numeric. The convergence threshold for the alignment objective
#'   function (Frobenius discrepancy).
#' @param add_names Logical. If `TRUE` (default), assigns Stan-format column names
#'   (e.g., `Lambda[1,1]`, `Lambda[2,1]`) to the output aligned matrix. If `FALSE`,
#'   the function checks whether `lambda_draws` already has column names: if it does,
#'   those names are preserved in the output; if it does not, `add_names` is
#'   automatically set to `TRUE`.
#'
#' @details
#' **Why Alignment is Necessary**
#'
#' Factor models are invariant to orthogonal transformations: if \eqn{\mathbf{\Lambda}}
#' is a valid loading matrix, then \eqn{\mathbf{\Lambda}^\bullet = \mathbf{\Lambda} \cdot \mathbf{Q}}
#' yields an identical model-implied covariance matrix for any orthogonal matrix
#' \eqn{\mathbf{Q}} where \eqn{\mathbf{Q} \cdot \mathbf{Q}^\top = \mathbf{I}_M}.
#' During MCMC sampling, this invariance creates a multimodal posterior where the
#' sampler explores equivalent, yet incompatible, rotational orientations.
#' Consequently, naive posterior summaries (such as averaging across draws) mix
#' these modes, yielding uninterpretable loading estimates that artificially cancel
#' out toward zero.
#'
#' **The Efficient RSP Algorithm (E-RSP)**
#'
#' This function implements the Efficient Rotation-Sign-Permutation (E-RSP)
#' algorithm (Rey-Sáez et al., 2026), an optimized version of the exact RSP
#' method originally proposed by Papastamoulis and Ntzoufras (2022).
#' While the original approach becomes computationally prohibitive as the
#' number of latent factors increases, E-RSP reduces the
#' alignment task to a Linear Assignment Problem (LAP). This
#' guarantees a globally optimal solution that scales to high-dimensional
#' models with negligible computational cost.
#'
#' The algorithm proceeds in two stages:
#'
#' 1. **Continuous Alignment (Varimax)**: Each posterior draw is rotated to a
#'    canonical simple structure to fix the continuous rotational degree of freedom.
#' 2. **Discrete Alignment (Signed-Permutation)**: The remaining discrete ambiguity
#'    (column permutations and sign flips) is resolved by mapping each draw to a
#'    common reference configuration.
#'
#' **Flattening Formats**
#'
#' The `lambda_draws` input must be a 2D matrix of dimensions S \eqn{\times} (J \eqn{\times} M),
#' where each row represents a flattened loading matrix from a single MCMC draw.
#' The `format` argument simply tells the function how your MCMC software flattened the
#' original \eqn{J \times M} loading matrix \eqn{\mathbf{\Lambda}} into that row vector.
#' The aligned output will be returned in this exact same flattened format.
#'
#' To illustrate, consider a \eqn{4 \times 3} loading matrix (4 items, 3 factors):
#' \deqn{\mathbf{\Lambda} = \begin{bmatrix}
#' \lambda_{11} & \lambda_{12} & \lambda_{13} \\
#' \lambda_{21} & \lambda_{22} & \lambda_{23} \\
#' \lambda_{31} & \lambda_{32} & \lambda_{33} \\
#' \lambda_{41} & \lambda_{42} & \lambda_{43}
#' \end{bmatrix}}
#'
#' Depending on the `format` specified, the input rows must be structured as follows:
#'
#' * **Column-major** (`format = "column_major"`): Elements are filled column by column
#'   (factor by factor). This corresponds to the `as.vector(Lambda)` behavior in R and
#'   is the default output layout for Stan:
#'   \deqn{\text{vec}(\mathbf{\Lambda}) = \big[ \underbrace{\lambda_{11}, \lambda_{21}, \lambda_{31}, \lambda_{41}}_{\text{Factor 1}}, \, \underbrace{\lambda_{12}, \lambda_{22}, \lambda_{32}, \lambda_{42}}_{\text{Factor 2}}, \, \underbrace{\lambda_{13}, \lambda_{23}, \lambda_{33}, \lambda_{43}}_{\text{Factor 3}} \big]}
#'
#' * **Row-major** (`format = "row_major"`): Elements are filled row by row (item by item).
#'   This corresponds to `as.vector(t(Lambda))` in R:
#'   \deqn{\text{vec}(\mathbf{\Lambda}^\top) = \big[ \underbrace{\lambda_{11}, \lambda_{12}, \lambda_{13}}_{\text{Item 1}}, \, \underbrace{\lambda_{21}, \lambda_{22}, \lambda_{23}}_{\text{Item 2}}, \, \underbrace{\lambda_{31}, \lambda_{32}, \lambda_{33}}_{\text{Item 3}}, \, \underbrace{\lambda_{41}, \lambda_{42}, \lambda_{43}}_{\text{Item 4}} \big]}
#'
#' @return A list containing:
#'   * `Lambda_hat_mcmc`: A numeric matrix of the aligned loading draws,
#'     maintaining the exact dimensions and format of the input `lambda_draws`.
#'     The attribute `nchains` is added to this matrix, indicating the number
#'     of chains (as passed in `n_chains`).
#'   * `Lambda_star`: A \eqn{J \times M} numeric matrix representing the final
#'     reference configuration (the posterior mean of the aligned draws).
#'   * `objective`: A numeric value indicating the final alignment objective
#'     (the total Frobenius discrepancy across all draws).
#'   * `sign_vectors`: A matrix of dimensions \eqn{S \times M} containing the
#'     sign flips applied to each draw.
#'   * `perm_vectors`: A matrix of dimensions \eqn{S \times M} containing the
#'     permutations applied to each draw.
#'   * `summary`: A `draws_summary` data frame computed via
#'     [posterior::summarise_draws()].
#'
#' @references
#' Papastamoulis, P., & Ntzoufras, I. (2022). On the identifiability of Bayesian
#' factor analytic models. *Statistics and Computing, 32*(2), 23.
#' <https://doi.org/10.1007/s11222-022-10084-4>
#'
#' Rey-Sáez, R., Vadillo, M. A., & Revuelta, J. (2026). *An Efficient
#' Rotation-Sign-Permutation Algorithm to Solve Rotational Indeterminacy in
#' Bayesian Exploratory Factor Analysis*. Manuscript submitted for publication.
#' <https://osf.io/5dutv/>
#'
#' @examples
#' \dontrun{
#' # Fit Bayesian EFA model
#' befa_fit <- befa(
#'   data = HS_data,
#'   n_factors = 3,
#'   rotate = "none",
#'   factor_scores = FALSE,
#'   compute_fit_indices = FALSE,
#'   compute_reliability = FALSE,
#'   backend = "cmdstanr",
#'   seed = 8
#' )
#'
#' # Extract unrotated posterior draws
#' lambda_unrotated <- extract_posterior_draws(befa_fit, pars = "Lambda")
#'
#' # See multimodality due to rotational indeterminacy
#' hist(lambda_unrotated[, 1],
#'   breaks = 100, col = "steelblue2",
#'   main = "Rotation indeterminacy", xlab = "Lambda[1,1]"
#' )
#'
#' # See that columns are ordered following a column-major order
#' # see Details on rsp_align function
#' colnames(lambda_unrotated)
#'
#' #  [1] "Lambda[1,1]" "Lambda[2,1]" "Lambda[3,1]" "Lambda[4,1]" "Lambda[5,1]" "Lambda[6,1]" "Lambda[7,1]"
#' #  [8] "Lambda[8,1]" "Lambda[9,1]" "Lambda[1,2]" "Lambda[2,2]" "Lambda[3,2]" "Lambda[4,2]" "Lambda[5,2]"
#' # [15] "Lambda[6,2]" "Lambda[7,2]" "Lambda[8,2]" "Lambda[9,2]" "Lambda[1,3]" "Lambda[2,3]" "Lambda[3,3]"
#' # [22] "Lambda[4,3]" "Lambda[5,3]" "Lambda[6,3]" "Lambda[7,3]" "Lambda[8,3]" "Lambda[9,3]"
#'
#' # Now, let's align posterior draws
#' lambda_aligned <- rsp_align(
#'   lambda_draws = lambda_unrotated,
#'   n_items = ncol(HS_data),
#'   n_factors = 3,
#'   n_chains = 4,
#'   format = "column_major"
#' )
#'
#' # Now, rotational indeterminacy has been solved!
#' hist(lambda_aligned$Lambda_hat_mcmc[, 1],
#'   breaks = 100, col = "steelblue2",
#'   main = "Aligned posterior distribution", xlab = "Lambda[1,1]"
#' )
#' }
#'
#' @export
rsp_align <- function(lambda_draws, n_items, n_factors,
                      n_chains = 1,
                      format = NULL,
                      max_iter = 1000, threshold = 1e-6,
                      add_names = TRUE) {
  # Validate format (required argument)
  if (is.null(format)) {
    stop(paste0(
      "Argument 'format' is required. ",
      "Please specify format = \"column_major\" or format = \"row_major\". ",
      "See the Details section of ?rsp_align for an illustration of both formats."
    ), call. = FALSE)
  }
  format <- match.arg(format, choices = c("column_major", "row_major"))

  # ────────────────────────── #
  #    Input control checks    #
  # ────────────────────────── #

  # Ensure consistent dimensions
  if (ncol(lambda_draws) != (n_items * n_factors)) {
    stop(sprintf(
      "Dimension mismatch: lambda_draws has %d columns, but n_items * n_factors = %d.",
      ncol(lambda_draws), n_items * n_factors
    ))
  }

  # Preserve existing column names if add_names = FALSE
  existing_names <- colnames(lambda_draws)
  if (!add_names) {
    if (is.null(existing_names)) {
      # No existing names: force add_names = TRUE
      add_names <- TRUE
    }
  }

  # From column-major to row-major
  if (format == "row_major") {
    # Column index order to transform Column-Major in Row-Major
    idx_row_to_col <- as.vector(matrix(seq_len(n_items * n_factors),
      nrow = n_items, ncol = n_factors,
      byrow = TRUE
    ))
    # Reorder draws
    lambda_draws <- lambda_draws[, idx_row_to_col]
  }

  # ────────────────────────────────────────────────── #
  #    Apply Orthogonal Rotation Criteria (per draw)   #
  # ────────────────────────────────────────────────── #

  if (n_factors > 1) {
    # Always use Varimax for alignment
    orth_l <- t(apply(lambda_draws, 1, function(x) {
      c(stats::varimax(x = matrix(x, ncol = n_factors), normalize = FALSE)$loadings)
    }))
  } else {
    # With M = 1, just use the original draws
    orth_l <- lambda_draws
  }

  # ─────────────────────────────────── #
  #    Define Reference Matrix          #
  # ─────────────────────────────────── #

  # Use the posterior mean as Papastamoulis & Ntzoufras (2022)
  lambda_ref <- matrix(colMeans(orth_l), nrow = n_items, ncol = n_factors)

  # ─────────────────────────────────── #
  #    Apply Efficient RSP algorithm    #
  # ─────────────────────────────────── #

  # Valid posterior orthogonal draws
  rsp_res <- rsp_exact_eff(
    Lambda_tilde = orth_l,
    Lambda_star_init = lambda_ref,
    J = n_items,
    M = n_factors,
    maxIter = max_iter,
    threshold = threshold
  )

  # ─────────────────────────────── #
  #    Force Positive Reflection    #
  # ─────────────────────────────── #

  # Calculate signs to force positive column means in Lambda_star
  s_final <- sign(colMeans(rsp_res$Lambda_star))

  # Update Lambda_star
  rsp_res$Lambda_star <- sweep(rsp_res$Lambda_star, 2, s_final, "*")

  # Update Lambda_hat_mcmc
  s_final_expanded <- rep(s_final, each = n_items)
  rsp_res$Lambda_hat_mcmc <- sweep(rsp_res$Lambda_hat_mcmc, 2, s_final_expanded, "*")

  # Update sign_vectors
  for (k in 1:n_factors) {
    if (s_final[k] == -1) {
      original_cols_for_k <- rsp_res$perm_vectors[, k]
      linear_indices <- (seq_len(nrow(rsp_res$sign_vectors))) + (original_cols_for_k - 1) * nrow(rsp_res$sign_vectors)
      rsp_res$sign_vectors[linear_indices] <- -1 * rsp_res$sign_vectors[linear_indices]
    }
  }

  if (format == "row_major") {
    # Reorder draws again to return in the same format
    rsp_res$Lambda_hat_mcmc <- rsp_res$Lambda_hat_mcmc[, order(idx_row_to_col)]
  }

  # Assign column names
  if (add_names) {
    # Assign Stan-format column names
    if (format == "column_major") {
      # Column-major: Lambda[1,1], Lambda[2,1], ..., Lambda[J,1], Lambda[1,2], ...
      col_names <- paste0("Lambda[", rep(seq_len(n_items), n_factors), ",", rep(seq_len(n_factors), each = n_items), "]")
    } else {
      # Row-major: Lambda[1,1], Lambda[1,2], ..., Lambda[1,M], Lambda[2,1], ...
      col_names <- paste0("Lambda[", rep(seq_len(n_items), each = n_factors), ",", rep(seq_len(n_factors), n_items), "]")
    }
    colnames(rsp_res$Lambda_hat_mcmc) <- col_names
  } else {
    # Preserve existing names from input
    colnames(rsp_res$Lambda_hat_mcmc) <- existing_names
  }

  # Add nchains attribute
  attr(rsp_res$Lambda_hat_mcmc, "nchains") <- n_chains

  # Finally, add posterior summaries
  rsp_res$summary <- posterior::summarise_draws(rsp_res$Lambda_hat_mcmc)

  return(rsp_res)
}
