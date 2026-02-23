#' Bayesian Reliability Estimation (McDonald's Omega)
#'
#' Computes model-based reliability coefficients (McDonald's Omega) from the
#' posterior distribution of a fitted BEFA model. By evaluating the reliability
#' analytically at each MCMC draw, it provides full posterior distributions
#' (and credible intervals) for both total-score and subscale reliabilities,
#' avoiding the need for asymptotic approximations.
#'
#' @param object A `befa` object returned by [befa()].
#' @param probs Numeric vector of length 2. Quantiles used to compute the
#'   credible intervals (default: `c(0.025, 0.975)` for a 95% interval).
#'
#' @details
#' ### Why McDonald's Omega?
#'
#' Unlike Cronbach's \eqn{\alpha}, which strictly assumes tau-equivalence
#' (equal factor loadings across all items), McDonald's \eqn{\omega} correctly
#' accounts for the congeneric factor structure estimated in Exploratory Factor
#' Analysis (EFA). Furthermore, by computing the coefficients draw-by-draw from
#' the posterior distribution, this function naturally propagates parameter
#' uncertainty, yielding full credible intervals for reliability estimates.
#'
#' ### Mathematical Formulation
#'
#' Let \eqn{J} be the number of items, \eqn{M} the number of factors,
#' \eqn{\mathbf{\Lambda}} the \eqn{J \times M} loading matrix, \eqn{\mathbf{\Psi}}
#' the diagonal matrix of unique variances (uniquenesses), and \eqn{\mathbf{\Sigma}}
#' the model-implied covariance matrix. \eqn{\mathbf{1}} is a vector of ones of length \eqn{J}.
#'
#' * **Omega Total (\eqn{\omega_t})**:
#'   Represents the proportion of variance in the unit-weighted total score
#'   (sum of all items) attributable to *all* common factors.
#'   \deqn{\large \omega_t = \Large{\frac{\sum_{m=1}^M (\sum_{j=1}^J \lambda_{jm})^2}{\mathbf{1}^\top \cdot \mathbf{\Sigma} \cdot \mathbf{1}}}}
#'
#' * **Omega Subscale (\eqn{\omega_{s,m}})**:
#'   Represents the reliability of the sum score when isolating the variance
#'   attributable to a specific factor \eqn{m}, treating the unique variances
#'   as measurement error.
#'   \deqn{\large \omega_{s,m} = \Large{\frac{(\sum_{j=1}^J \lambda_{jm})^2}{(\sum_{j=1}^J \lambda_{jm})^2 + \sum_{j=1}^J \psi_j}}}
#'
#' @return A list of class `befa_reliability` containing:
#'   * `omega_draws`: A `draws_array` (iterations \eqn{\times} chains \eqn{\times} P)
#'     containing the posterior draws of all reliability coefficients (`omega_total`,
#'     `omega_F1`, ..., `omega_FM`).
#'   * `summary`: A `draws_summary` data frame computed via
#'     [posterior::summarise_draws()].
#'
#' @references
#' McDonald, R. P. (1999). *Test Theory: A Unified Treatment*. Lawrence Erlbaum Associates.
#'
#' Zinbarg, R. E., Revelle, W., Yovel, I., & Li, W. (2005). Cronbach's \eqn{\alpha},
#' Revelle's \eqn{\beta}, and McDonald's \eqn{\omega_H}: Their relations with
#' each other and two alternative conceptualizations of reliability.
#' *Psychometrika, 70*(1), 123-133. <https://doi.org/10.1007/s11336-003-0974-7>
#'
#' @examples
#' \dontrun{
#' # ──────────────────────────────────────────────────────────────────────── #
#' #    1. Bayesian Reliability estimates after fitting Bayesian EFA model    #
#' # ──────────────────────────────────────────────────────────────────────── #
#'
#' # Fit Bayesian EFA model
#' befa_fit <- befa(
#'   data = HS_data,
#'   n_factors = 3,
#'   factor_scores = FALSE,
#'   compute_fit_indices = FALSE,
#'   compute_reliability = FALSE,
#'   backend = "cmdstanr",
#'   seed = 2026
#' )
#'
#' # Compute Bayesian Reliability estimates from a fitted befa object
#' bayesian_omega <- befa_reliability(befa_fit, probs = c(.025, .975))
#'
#' # Posterior summaries
#' bayesian_omega$summary
#'
#' # # A tibble: 4 × 10
#' #   variable     mean median     sd    mad    q5   q95  rhat ess_bulk ess_tail
#' #   <chr>       <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
#' # 1 omega_total 0.843  0.844 0.0108 0.0106 0.825 0.860 1.00     4628.    3442.
#' # 2 omega_F1    0.545  0.550 0.0449 0.0422 0.465 0.609 1.000    3702.    3712.
#' # 3 omega_F2    0.727  0.728 0.0187 0.0172 0.697 0.755 1.00     4126.    3109.
#' # 4 omega_F3    0.599  0.601 0.0387 0.0381 0.534 0.659 1.00     4248.    3413.
#'
#' # Posterior draws (transform to matrix)
#' posterior::as_draws_matrix(bayesian_omega$omega_draws)
#'
#' # # A draws_matrix: 1000 iterations, 4 chains, and 4 variables
#' #     variable
#' # draw omega_total omega_F1 omega_F2 omega_F3
#' #   1         0.84     0.43     0.75     0.63
#' #   2         0.82     0.51     0.69     0.59
#' #   3         0.85     0.54     0.72     0.64
#' #   4         0.85     0.55     0.72     0.65
#' #   5         0.83     0.51     0.72     0.57
#' #   6         0.83     0.58     0.70     0.53
#' #   7         0.84     0.56     0.70     0.59
#' #   8         0.84     0.52     0.73     0.61
#' #   9         0.86     0.55     0.75     0.64
#' #   10        0.85     0.62     0.73     0.59
#' # # ... with 3990 more draws
#'
#' # We can plot the posterior distribution of reliability estimates
#' par(mfrow = c(2, 2))
#' hist(posterior::as_draws_matrix(bayesian_omega$omega_draws)[, "omega_total"],
#'   breaks = 100, col = "grey90",
#'   main = "Histogram of omega index: Full scale",
#'   xlab = "Omega"
#' )
#' abline(
#'   v = bayesian_omega$summary$mean[1], col = "firebrick2",
#'   lwd = 3, lty = 2
#' )
#' hist(posterior::as_draws_matrix(bayesian_omega$omega_draws)[, "omega_F1"],
#'   breaks = 100, col = "grey90",
#'   main = "Histogram of omega index: Factor 1",
#'   xlab = "Omega"
#' )
#' abline(
#'   v = bayesian_omega$summary$mean[2], col = "firebrick2",
#'   lwd = 3, lty = 2
#' )
#' hist(posterior::as_draws_matrix(bayesian_omega$omega_draws)[, "omega_F2"],
#'   breaks = 100, col = "grey90",
#'   main = "Histogram of omega index: Factor 2",
#'   xlab = "Omega"
#' )
#' abline(
#'   v = bayesian_omega$summary$mean[3], col = "firebrick2",
#'   lwd = 3, lty = 2
#' )
#' hist(posterior::as_draws_matrix(bayesian_omega$omega_draws)[, "omega_F3"],
#'   breaks = 100, col = "grey90",
#'   main = "Histogram of omega index: Factor 3",
#'   xlab = "Omega"
#' )
#' abline(
#'   v = bayesian_omega$summary$mean[4], col = "firebrick2",
#'   lwd = 3, lty = 2
#' )
#' par(mfrow = c(1, 1))
#'
#' # ─────────────────────────────────────────────────────────────────────────────
#'
#' # ─────────────────────────────────────────────────────────────────────── #
#' #    2. Bayesian Reliability estimates when fitting Bayesian EFA model    #
#' # ─────────────────────────────────────────────────────────────────────── #
#'
#' # Fit Bayesian model with compute_fit_indices = TRUE
#' befa_fit <- befa(
#'   data = HS_data,
#'   n_factors = 3,
#'   factor_scores = FALSE,
#'   compute_fit_indices = FALSE,
#'   compute_reliability = TRUE,
#'   backend = "cmdstanr",
#'   seed = 8
#' )
#'
#' # befa_reliability output is inside the "reliability" object
#' # befa_fit$reliability
#'
#' # ─────────────────────────────────────────────────────────────────────────────
#' }
#'
#' @export
befa_reliability <- function(object, probs = c(0.025, 0.975)) {
  # Return cached results if already computed
  if (!is.null(object$reliability)) {
    return(object$reliability)
  }

  J <- object$stan_data$J
  M <- object$n_factors
  # Extract posterior draws preserving chain structure
  post_draws <- posterior::as_draws_array(object$stanfit)
  n_iter <- dim(post_draws)[1]
  n_chains <- dim(post_draws)[2]
  var_names <- dimnames(post_draws)[[3]]

  # Lambda and Psi are always on the same scale (correlation for UV/uni,
  # raw for normal). Because omega is scale-invariant, we compute
  # Sigma_implied = Lambda Lambda' + diag(Psi) directly.
  lambda_idx <- grep("^Lambda", var_names)
  psi_idx <- grep("^Psi", var_names)

  # Parameter names for the omega array
  omega_names <- c("omega_total", paste0("omega_F", seq_len(M)))
  n_params <- length(omega_names)

  # Output array: iterations x chains x (M+1)
  omega_array <- array(NA_real_,
    dim = c(n_iter, n_chains, n_params),
    dimnames = list(
      iteration = seq_len(n_iter),
      chain = seq_len(n_chains),
      variable = omega_names
    )
  )

  # Compute reliability measures draw-by-draw, preserving chain structure
  for (ch in seq_len(n_chains)) {
    for (s in seq_len(n_iter)) {
      L_s <- matrix(as.numeric(post_draws[s, ch, lambda_idx]), nrow = J, ncol = M)
      psi_s <- as.numeric(post_draws[s, ch, psi_idx])

      # Model-implied covariance: Sigma = LL' + diag(Psi)
      Sigma_implied <- tcrossprod(L_s) + diag(psi_s)

      # Omega per factor (i.e., subscale reliability)
      for (m in seq_len(M)) {
        num_m <- sum(L_s[, m])^2
        den_m <- num_m + sum(psi_s)
        omega_array[s, ch, m + 1] <- num_m / den_m
      }

      # Omega total (i.e., full scale reliability)
      common_var_total <- sum(colSums(L_s)^2)
      omega_array[s, ch, 1] <- common_var_total / sum(Sigma_implied)
    }
  }

  # Convert to posterior draws_array
  omega_draws <- posterior::as_draws_array(omega_array)

  # Compute summary with credible intervals based on probs argument
  summ <- posterior::summarise_draws(omega_draws)

  # ── Build convenience fields (total / subscales) ──
  # Compute quantiles per variable using probs
  omega_mat <- posterior::as_draws_matrix(omega_draws)

  # Total reliability (omega_total)
  total_draws <- as.numeric(omega_mat[, "omega_total"])
  qs_total <- stats::quantile(total_draws, probs = probs, na.rm = TRUE)
  total <- list(
    Estimate  = mean(total_draws, na.rm = TRUE),
    SD        = stats::sd(total_draws, na.rm = TRUE),
    Lower_95  = unname(qs_total[1]),
    Upper_95  = unname(qs_total[2])
  )

  # Subscale reliability (omega_F1, omega_F2, ...)
  sub_names <- paste0("omega_F", 1:M)
  sub_est <- numeric(M)
  sub_sd <- numeric(M)
  sub_lo <- numeric(M)
  sub_hi <- numeric(M)
  for (m in seq_len(M)) {
    sub_draws <- as.numeric(omega_mat[, sub_names[m]])
    qs_sub <- stats::quantile(sub_draws, probs = probs, na.rm = TRUE)
    sub_est[m] <- mean(sub_draws, na.rm = TRUE)
    sub_sd[m] <- stats::sd(sub_draws, na.rm = TRUE)
    sub_lo[m] <- unname(qs_sub[1])
    sub_hi[m] <- unname(qs_sub[2])
  }
  subscales <- data.frame(
    Estimate  = sub_est,
    SD        = sub_sd,
    Lower_95  = sub_lo,
    Upper_95  = sub_hi
  )
  rownames(subscales) <- paste0("F", 1:M)

  # Assembly of final results
  res <- list(
    omega_draws = omega_draws,
    summary     = summ,
    total       = total,
    subscales   = subscales
  )

  class(res) <- "befa_reliability"
  return(res)
}
