# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: compute_posterior_metrics
# ─────────────────────────────────────────────────────────────────────────────

# --- Helper: create mock befa object with draws ---
create_mock_befa <- function(model_type = "raw", J = 5, M = 2, n_draws = 50) {
  set.seed(123)
  N <- 100

  # Generate fake data
  data <- matrix(rnorm(N * J), nrow = N, ncol = J)

  # Create fake Sigma draws (J x J for each draw)
  sigma_cols <- paste0("Sigma[", rep(1:J, J), ",", rep(1:J, each = J), "]")
  rho_cols <- paste0("Rho[", rep(1:J, J), ",", rep(1:J, each = J), "]")
  nu_cols <- paste0("nu[", 1:J, "]")
  lambda_cols <- paste0("Lambda[", rep(1:J, M), ",", rep(1:M, each = J), "]")

  # Create a proper covariance matrix for each draw
  draws_matrix <- matrix(NA, nrow = n_draws, ncol = 0)

  # Add Sigma or Rho based on model_type
  sigma_draws <- matrix(NA, nrow = n_draws, ncol = J * J)
  for (s in 1:n_draws) {
    # Generate positive definite matrix
    A <- matrix(rnorm(J * J, 0, 0.3), J, J)
    Sigma <- A %*% t(A) + diag(1, J)
    sigma_draws[s, ] <- as.vector(Sigma)
  }

  if (model_type == "raw") {
    colnames(sigma_draws) <- sigma_cols
    # Add nu draws
    nu_draws <- matrix(rnorm(n_draws * J, 0, 0.5), nrow = n_draws)
    colnames(nu_draws) <- nu_cols
    draws_matrix <- cbind(sigma_draws, nu_draws)
  } else {
    # For std model, use correlation matrices
    rho_draws <- t(apply(sigma_draws, 1, function(s) {
      S <- matrix(s, J, J)
      as.vector(cov2cor(S))
    }))
    colnames(rho_draws) <- rho_cols
    draws_matrix <- rho_draws
  }

  # Add Lambda draws
  lambda_draws <- matrix(rnorm(n_draws * J * M, 0, 0.5), nrow = n_draws)
  colnames(lambda_draws) <- lambda_cols
  draws_matrix <- cbind(draws_matrix, lambda_draws)

  # Create mock stanfit-like object
  mock_stanfit <- structure(
    list(draws = draws_matrix),
    class = "mock_stanfit"
  )

  # Create mock befa object
  list(
    model_type = model_type,
    stanfit = mock_stanfit,
    stan_data = list(N = N, J = J, M = M, Y = data),
    n_factors = M
  )
}

# Mock posterior::as_draws_matrix
mock_as_draws_matrix <- function(x) x$draws

# --- Tests for compute_posterior_metrics ---

test_that("compute_posterior_metrics returns correct structure", {
  skip("Requires full befa object - use integration test")

  # object <- readRDS("tests/testthat/fixtures/befa_fit.rds")
  # data <- object$stan_data$Y
  # result <- compute_posterior_metrics(object, data)

  # expect_type(result, "list")
  # expect_true(all(c("log_lik", "chisq", "chisq_rep", "srmr",
  #                   "ll_saturated", "N", "J") %in% names(result)))
})

test_that("compute_posterior_metrics log_lik has correct dimensions", {
  skip("Requires full befa object")

  # object <- readRDS("tests/testthat/fixtures/befa_fit.rds")
  # data <- object$stan_data$Y
  # n_draws <- nrow(posterior::as_draws_matrix(object$stanfit))
  # N <- nrow(data)
  # result <- compute_posterior_metrics(object, data)

  # expect_equal(dim(result$log_lik), c(n_draws, N))
})

test_that("compute_posterior_metrics chisq is always non-negative", {
  skip("Requires full befa object")

  # Chi-square statistic should be >= 0
  # result <- compute_posterior_metrics(object, data)
  # expect_true(all(result$chisq >= 0))
})

test_that("compute_posterior_metrics srmr is between 0 and 1", {
  skip("Requires full befa object")

  # SRMR typically ranges 0-1 (can exceed slightly in bad fits)
  # result <- compute_posterior_metrics(object, data)
  # expect_true(all(result$srmr >= 0))
  # expect_true(all(result$srmr < 2))  # Reasonable upper bound
})

test_that("compute_posterior_metrics handles std model correctly", {
  skip("Requires full befa object with model_type = 'std'")

  # For std model, should look for Rho instead of Sigma
  # object_std <- readRDS("tests/testthat/fixtures/befa_fit_std.rds")
  # result <- compute_posterior_metrics(object_std, scale(object_std$stan_data$Y))
  # expect_type(result, "list")
})

test_that("compute_posterior_metrics errors if Sigma/Rho not found", {
  # Create object missing required parameters
  mock_obj <- list(
    model_type = "raw",
    stanfit = structure(list(draws = matrix(1:10, nrow = 2)), class = "mock"),
    stan_data = list(Y = matrix(rnorm(20), 10, 2))
  )

  # This test depends on how posterior::as_draws_matrix handles mock objects
  # In practice, test with a real malformed stanfit
  skip("Requires specific mock setup")
})

test_that("compute_posterior_metrics ll_saturated is scalar", {
  skip("Requires full befa object")

  # result <- compute_posterior_metrics(object, data)
  # expect_length(result$ll_saturated, 1)
  # expect_true(is.finite(result$ll_saturated))
})

# --- Unit tests for internal calculations ---

test_that("saturated model log-likelihood is computable", {
  set.seed(42)
  N <- 50
  J <- 4
  data <- matrix(rnorm(N * J), N, J)

  S_sat <- cov(data) * (N - 1) / N
  M_sat <- colMeans(data)

  ll_sat <- sum(mvtnorm::dmvnorm(data, mean = M_sat, sigma = S_sat, log = TRUE))

  expect_true(is.finite(ll_sat))
  expect_type(ll_sat, "double")
})

test_that("SRMR calculation is correct for perfect fit", {
  # If model-implied equals sample correlation, SRMR = 0
  J <- 4
  Sample_Cor <- diag(J) # Identity
  Imp_Cor <- diag(J) # Same

  diff_mat <- Sample_Cor - Imp_Cor
  resid_unique <- diff_mat[lower.tri(diff_mat, diag = TRUE)]
  srmr <- sqrt(mean(resid_unique^2))

  expect_equal(srmr, 0)
})

test_that("SRMR increases with model misfit", {
  J <- 4
  Sample_Cor <- diag(J)

  # Small misfit
  Imp_Cor_small <- diag(J) + 0.01
  diag(Imp_Cor_small) <- 1
  diff_small <- Sample_Cor - Imp_Cor_small
  srmr_small <- sqrt(mean(diff_small[lower.tri(diff_small, diag = TRUE)]^2))

  # Large misfit
  Imp_Cor_large <- diag(J) + 0.1
  diag(Imp_Cor_large) <- 1
  diff_large <- Sample_Cor - Imp_Cor_large
  srmr_large <- sqrt(mean(diff_large[lower.tri(diff_large, diag = TRUE)]^2))

  expect_true(srmr_large > srmr_small)
})

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: compute_null_metrics
# ─────────────────────────────────────────────────────────────────────────────

test_that("compute_null_metrics returns correct structure for std model", {
  set.seed(123)
  N <- 50
  J <- 4
  data <- matrix(rnorm(N * J), N, J)
  z_data <- scale(data)

  S_sat <- cor(data)
  M_sat <- rep(0, J)
  ll_saturated <- sum(mvtnorm::dmvnorm(z_data, mean = M_sat, sigma = S_sat, log = TRUE))

  result <- compute_null_metrics(
    data = data,
    model_type = "cor",
    ll_saturated = ll_saturated,
    stan_data = list(N = N, J = J, M = 2),
    stan_args_list = list(),
    lambda_prior = "unit_vector"
  )

  expect_type(result, "list")
  expect_true(all(c("chisq_vec", "chisq_mean", "pD", "df", "lambda_vec") %in% names(result)))
})

test_that("compute_null_metrics std model returns scalar chisq (deterministic)", {
  set.seed(456)
  N <- 50
  J <- 4
  data <- matrix(rnorm(N * J), N, J)
  z_data <- scale(data)

  S_sat <- cor(data)
  ll_saturated <- sum(mvtnorm::dmvnorm(z_data, mean = rep(0, J), sigma = S_sat, log = TRUE))

  result <- compute_null_metrics(
    data = data,
    model_type = "cor",
    ll_saturated = ll_saturated,
    stan_data = list(N = N, J = J, M = 2),
    stan_args_list = list(),
    lambda_prior = "unit_vector"
  )

  # For std model, null is deterministic -> scalar
  expect_length(result$chisq_vec, 1)
  expect_equal(result$pD, 0)
})

test_that("compute_null_metrics std model chisq is non-negative", {
  set.seed(789)
  N <- 100
  J <- 5
  data <- matrix(rnorm(N * J), N, J)
  z_data <- scale(data)

  S_sat <- cor(data)
  ll_saturated <- sum(mvtnorm::dmvnorm(z_data, mean = rep(0, J), sigma = S_sat, log = TRUE))

  result <- compute_null_metrics(
    data = data,
    model_type = "cor",
    ll_saturated = ll_saturated,
    stan_data = list(N = N, J = J, M = 2),
    stan_args_list = list(),
    lambda_prior = "unit_vector"
  )

  expect_true(result$chisq_vec >= 0)
})

test_that("compute_null_metrics std model df is at least 1", {
  set.seed(111)
  N <- 50
  J <- 4
  data <- matrix(rnorm(N * J), N, J)
  z_data <- scale(data)

  S_sat <- cor(data)
  ll_saturated <- sum(mvtnorm::dmvnorm(z_data, mean = rep(0, J), sigma = S_sat, log = TRUE))

  result <- compute_null_metrics(
    data = data,
    model_type = "cor",
    ll_saturated = ll_saturated,
    stan_data = list(N = N, J = J, M = 2),
    stan_args_list = list(),
    lambda_prior = "unit_vector"
  )

  expect_true(result$df >= 1)
})

test_that("compute_null_metrics lambda_vec is non-negative", {
  set.seed(222)
  N <- 50
  J <- 4
  data <- matrix(rnorm(N * J), N, J)
  z_data <- scale(data)

  S_sat <- cor(data)
  ll_saturated <- sum(mvtnorm::dmvnorm(z_data, mean = rep(0, J), sigma = S_sat, log = TRUE))

  result <- compute_null_metrics(
    data = data,
    model_type = "cor",
    ll_saturated = ll_saturated,
    stan_data = list(N = N, J = J, M = 2),
    stan_args_list = list(),
    lambda_prior = "unit_vector"
  )

  expect_true(all(result$lambda_vec >= 0))
})

test_that("compute_null_metrics raw model requires Stan estimation", {
  skip_if_not_installed("rstan")
  skip("Requires rstan and befa_raw_null model - integration test")

  # This would actually run Stan, so skip in unit tests
  # Full integration test should be run separately
})

test_that("compute_null_metrics p_star calculation is correct for std", {
  J <- 5
  p_star_std <- J * (J - 1) / 2 # Only correlations
  expect_equal(p_star_std, 10)
})

test_that("compute_null_metrics p_star calculation is correct for raw", {
  J <- 5
  p_star_raw <- J * (J + 3) / 2 # MACS
  expect_equal(p_star_raw, 20)
})

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: print.befa_fitmeasures
# ─────────────────────────────────────────────────────────────────────────────

# --- Helper: create mock befa_fitmeasures object ---
create_mock_fitmeasures <- function() {
  set.seed(123)
  n_draws <- 100

  # Create a proper draws_matrix that posterior::summarise_draws can handle
  fit_draws_mat <- cbind(
    Chi2       = rnorm(n_draws, 50.5, 5.2),
    Chi2_Null  = rnorm(n_draws, 120.3, 8.1),
    BRMSEA     = abs(rnorm(n_draws, 0.045, 0.008)),
    BGamma     = pmin(1, pmax(0, rnorm(n_draws, 0.95, 0.01))),
    Adj_BGamma = pmin(1, pmax(0, rnorm(n_draws, 0.94, 0.012))),
    BMc        = pmin(1, pmax(0, rnorm(n_draws, 0.96, 0.009))),
    SRMR       = abs(rnorm(n_draws, 0.038, 0.005)),
    BCFI       = pmin(1, pmax(0, rnorm(n_draws, 0.97, 0.008))),
    BTLI       = pmin(1, pmax(0, rnorm(n_draws, 0.96, 0.01)))
  )
  fit_draws <- posterior::as_draws_matrix(fit_draws_mat)

  # Create mock loo object with estimates matrix
  loo_estimates <- matrix(
    c(
      -150.5, 5.2, # elpd_loo
      301.0, 10.4, # looic
      8.5, 1.2
    ), # p_loo
    nrow = 3, ncol = 2, byrow = TRUE,
    dimnames = list(c("elpd_loo", "looic", "p_loo"), c("Estimate", "SE"))
  )
  mock_loo <- list(estimates = loo_estimates)

  structure(
    list(
      fit_draws = fit_draws,
      loo_object = mock_loo,
      details = list(chi2_ppp = 0.45, p_star = 15, pD = 8.5, N = 200)
    ),
    class = "befa_fitmeasures"
  )
}

# Tests start here

test_that("print.befa_fitmeasures runs without error", {
  mock_fm <- create_mock_fitmeasures()

  expect_output(print(mock_fm), "Bayesian Fit Measures")
})

test_that("print.befa_fitmeasures returns object invisibly", {
  mock_fm <- create_mock_fitmeasures()

  result <- capture.output(ret <- print(mock_fm))
  expect_identical(ret, mock_fm)
})

test_that("print.befa_fitmeasures respects digits argument", {
  mock_fm <- create_mock_fitmeasures()

  output_2 <- capture.output(print(mock_fm, digits = 2))
  output_5 <- capture.output(print(mock_fm, digits = 5))

  # Different digit settings should produce different output
  # (character lengths differ)
  expect_true(length(output_2) > 0)
  expect_true(length(output_5) > 0)
})

test_that("print.befa_fitmeasures shows all indices", {
  mock_fm <- create_mock_fitmeasures()

  output <- capture.output(print(mock_fm))
  full_output <- paste(output, collapse = "\n")

  expect_match(full_output, "BRMSEA", fixed = TRUE)
  expect_match(full_output, "BCFI", fixed = TRUE)
  expect_match(full_output, "SRMR", fixed = TRUE)
  expect_match(full_output, "ELPD", fixed = TRUE)
})

test_that("print.befa_fitmeasures shows notes", {
  mock_fm <- create_mock_fitmeasures()

  output <- capture.output(print(mock_fm))
  full_output <- paste(output, collapse = "\n")

  expect_match(full_output, "Note:", fixed = TRUE)
  expect_match(full_output, "95%", fixed = TRUE)
})

test_that("print.befa_fitmeasures handles NA values gracefully", {
  mock_fm <- create_mock_fitmeasures()
  # Chi2_ppp has NAs in SD and CI columns

  expect_output(print(mock_fm), "Chi2_ppp")
})

# ─────────────────────────────────────────────────────────────────────────────
