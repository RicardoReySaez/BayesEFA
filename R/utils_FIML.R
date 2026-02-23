#' Compute FIML Log-Likelihood (R equivalent of Stan's incomplete_obs_lpdf)
#'
#' For each observation, computes the MVN log-density using only observed
#' variables (marginalizing over missing ones). This is the proper FIML approach.
#' It is highly optimized to compute complete cases vectorized.
#'
#' @param data_fiml A list containing `Y_comp` (matrix of complete cases),
#'   `Y_miss` (matrix of incomplete cases), `comp_idx` (original indices for complete cases),
#'   and `miss_idx` (original indices for incomplete cases).
#' @param mu Vector. Mean vector (length J).
#' @param Sigma Matrix. Covariance matrix (J x J).
#' @return Numeric vector of log-likelihoods matching the original N.
#' @keywords internal
#' @noRd
compute_fiml_loglik <- function(data_fiml, mu, Sigma) {
  # Total number of observations is the sum of complete and incomplete cases
  N_comp <- length(data_fiml$comp_idx)
  N_miss <- length(data_fiml$miss_idx)
  N_total <- N_comp + N_miss

  ll_vec <- numeric(N_total)

  # 1. Fast-Path for complete cases (Vectorized C++ call)
  if (N_comp > 0) {
    ll_vec[data_fiml$comp_idx] <- mvnfast::dmvn(data_fiml$Y_comp, mu = mu, sigma = Sigma, log = TRUE)
  }

  # 2. Row-by-row FIML for incomplete cases
  if (N_miss > 0) {
    Y_miss <- data_fiml$Y_miss

    # Pre-allocate subset vector for missing cases
    ll_miss <- numeric(N_miss)

    for (i in seq_len(N_miss)) {
      obs_idx <- which(!is.na(Y_miss[i, ]))
      n_obs <- length(obs_idx)

      if (n_obs == 0) {
        ll_miss[i] <- 0 # No observed data -> no contribution
      } else {
        # Marginalize over observed dimensions
        y_obs <- Y_miss[i, obs_idx]
        mu_obs <- mu[obs_idx]
        Sigma_obs <- Sigma[obs_idx, obs_idx, drop = FALSE]

        ll_miss[i] <- mvnfast::dmvn(matrix(y_obs, nrow = 1), mu = mu_obs, sigma = Sigma_obs, log = TRUE)
      }
    }
    # Assign the evaluated missing case likelihoods back to their original indices
    ll_vec[data_fiml$miss_idx] <- ll_miss
  }

  return(ll_vec)
}


# ═════════════════════════════════════════════════════════════════════════════════════
# EM Algorithm for Saturated Model Moments with Missing Data
# Based on lavaan's lav_mvnorm_missing_h1_estimate_moments and lav_mvnorm_missing_estep
# See: lavaan/R/lav_mvnorm_missing_h1.R by Yves Rosseel
# ═════════════════════════════════════════════════════════════════════════════════════

#' E-Step for EM Algorithm (Internal)
#'
#' Computes expected sufficient statistics T1 (sum) and T2 (sum of squares)
#' for data with missing values.
#'
#' @param Y Matrix with NAs.
#' @param Mu Current mean estimate.
#' @param Sigma Current covariance estimate.
#' @return List with T1 (vector) and T2 (matrix).
#' @keywords internal
#' @noRd
em_E_step <- function(Y, Mu, Sigma) {
  N <- nrow(Y)
  J <- ncol(Y)
  T1 <- numeric(J)
  T2 <- matrix(0, J, J)

  for (i in seq_len(N)) {
    obs_idx <- which(!is.na(Y[i, ]))
    mis_idx <- which(is.na(Y[i, ]))

    if (length(mis_idx) == 0) {
      # Complete case: direct contribution
      T1 <- T1 + Y[i, ]
      T2 <- T2 + tcrossprod(Y[i, ])
    } else if (length(obs_idx) > 0) {
      # Incomplete case: impute missing values
      Sigma_22_inv <- solve(Sigma[obs_idx, obs_idx, drop = FALSE])
      Sigma_12 <- Sigma[mis_idx, obs_idx, drop = FALSE]

      # Conditional mean: E[Y_miss | Y_obs]
      y_obs <- Y[i, obs_idx]
      y_imp <- Mu[mis_idx] + Sigma_12 %*% Sigma_22_inv %*% (y_obs - Mu[obs_idx])

      # Complete the observation
      y_complete <- numeric(J)
      y_complete[obs_idx] <- y_obs
      y_complete[mis_idx] <- as.vector(y_imp)

      T1 <- T1 + y_complete

      # T2 with conditional covariance correction
      T2_i <- tcrossprod(y_complete)
      # Correction: Var[Y_miss | Y_obs] = Sigma_11 - Sigma_12 * Sigma_22^{-1} * Sigma_21
      cond_cov <- Sigma[mis_idx, mis_idx, drop = FALSE] - Sigma_12 %*% Sigma_22_inv %*% t(Sigma_12)
      T2_i[mis_idx, mis_idx] <- T2_i[mis_idx, mis_idx] + cond_cov
      T2 <- T2 + T2_i
    }
  }

  list(T1 = T1, T2 = T2)
}


#' Estimate Saturated Model Moments via EM (Internal)
#'
#' Estimates mean vector and covariance matrix from data with missing values
#' using the EM algorithm. Ensures valid chi-square comparisons with FIML.
#'
#' @param Y Matrix with NAs.
#' @param max_iter Maximum EM iterations.
#' @param tol Convergence tolerance.
#' @return List with Mu and Sigma.
#' @keywords internal
#' @noRd
estimate_h1_moments_em <- function(Y, max_iter = 500L, tol = 1e-05) {
  N <- nrow(Y)

  # Initialize: column means (ignoring NAs), diagonal covariance
  Mu <- colMeans(Y, na.rm = TRUE)
  Yc <- sweep(Y, 2, Mu)
  var0 <- colMeans(Yc^2, na.rm = TRUE)
  var0[!is.finite(var0) | var0 == 0] <- 1
  Sigma <- diag(var0)

  for (iter in seq_len(max_iter)) {
    Mu0 <- Mu
    Sigma0 <- Sigma

    # E-step
    E_step <- em_E_step(Y, Mu, Sigma)

    # M-step
    Mu <- E_step$T1 / N
    Sigma <- E_step$T2 / N - tcrossprod(Mu)

    # Ensure positive definite (ridge if needed)
    ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev < 1e-6)) {
      diag(Sigma) <- diag(Sigma) + max(diag(Sigma)) * 1e-08
    }

    # Convergence check
    delta <- max(abs(c(Mu - Mu0, as.vector(Sigma - Sigma0))))
    if (delta < tol) break
  }

  list(Mu = Mu, Sigma = Sigma)
}
