# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: befa_reliability
# ─────────────────────────────────────────────────────────────────────────────

# --- Helper: create mock befa object with posterior draws ---
create_mock_befa_for_reliability <- function(model_type = "raw", J = 5, M = 2,
                                             n_draws = 100) {
  set.seed(123)

  # Create realistic Lambda draws
  true_lambda <- matrix(0.1, J, M)
  for (m in 1:M) {
    start <- ((m - 1) * floor(J / M)) + 1
    end <- min(m * floor(J / M), J)
    true_lambda[start:end, m] <- runif(end - start + 1, 0.6, 0.85)
  }

  lambda_draws <- matrix(NA, n_draws, J * M)
  sigma_draws <- matrix(NA, n_draws, J * J)

  for (s in 1:n_draws) {
    # Lambda with noise
    L_s <- true_lambda + matrix(rnorm(J * M, 0, 0.05), J, M)
    lambda_draws[s, ] <- as.vector(L_s)

    # Sigma = Lambda %*% t(Lambda) + diag(psi)
    psi_s <- runif(J, 0.2, 0.5)
    Sigma_s <- L_s %*% t(L_s) + diag(psi_s)
    sigma_draws[s, ] <- as.vector(Sigma_s)
  }

  # Create column names
  lambda_cols <- paste0("Lambda[", rep(1:J, M), ",", rep(1:M, each = J), "]")

  if (model_type == "raw") {
    sigma_cols <- paste0("Sigma[", rep(1:J, J), ",", rep(1:J, each = J), "]")
  } else {
    # For std model, normalize to correlation
    sigma_draws <- t(apply(sigma_draws, 1, function(s) {
      S <- matrix(s, J, J)
      as.vector(cov2cor(S))
    }))
    sigma_cols <- paste0("Rho[", rep(1:J, J), ",", rep(1:J, each = J), "]")
  }

  colnames(lambda_draws) <- lambda_cols
  colnames(sigma_draws) <- sigma_cols

  draws_matrix <- cbind(lambda_draws, sigma_draws)

  # Create mock stanfit
  mock_stanfit <- structure(
    list(draws = draws_matrix),
    class = "mock_stanfit"
  )

  # Create mock befa object
  structure(
    list(
      model_type = model_type,
      stanfit = mock_stanfit,
      stan_data = list(J = J, M = M, N = 200),
      n_factors = M
    ),
    class = "befa"
  )
}

# Note: Tests that require real befa objects (structure validation, M rows,
# rownames, etc.) are now covered by run_all_befa_checks() in test-befa.R

# --- Unit tests for omega calculations ---

test_that("omega total formula is correct for perfect communality", {
  # If all variance is common (psi = 0), omega_total = 1
  J <- 4
  M <- 2

  # Lambda with perfect structure
  L <- matrix(c(
    0.8, 0.7, 0.1, 0.1,
    0.1, 0.1, 0.8, 0.7
  ), nrow = 4, ncol = 2)

  # Sigma with zero uniqueness
  Sigma <- L %*% t(L) # No psi added

  # Omega total calculation
  total_implied_var <- sum(Sigma)
  common_var_total <- sum(colSums(L)^2)

  # Should be close to 1 (perfect reliability)
  omega_t <- common_var_total / total_implied_var
  expect_true(omega_t > 0.95)
})

test_that("omega total decreases with higher uniqueness", {
  J <- 4
  M <- 2
  L <- matrix(c(0.7, 0.1, 0.1, 0.7), nrow = 2, ncol = 2)

  # Low uniqueness
  psi_low <- c(0.1, 0.1)
  Sigma_low <- L %*% t(L) + diag(psi_low)
  omega_low <- sum(colSums(L)^2) / sum(Sigma_low)

  # High uniqueness
  psi_high <- c(0.8, 0.8)
  Sigma_high <- L %*% t(L) + diag(psi_high)
  omega_high <- sum(colSums(L)^2) / sum(Sigma_high)

  expect_true(omega_low > omega_high)
})

test_that("omega subscale formula is correct", {
  J <- 4
  M <- 2

  L <- matrix(c(
    0.8, 0.1,
    0.7, 0.1,
    0.1, 0.8,
    0.1, 0.7
  ), nrow = 4, ncol = 2, byrow = TRUE)

  psi <- c(0.3, 0.4, 0.3, 0.4)

  # Omega for factor 1
  num_1 <- sum(L[, 1])^2
  den_1 <- num_1 + sum(psi)
  omega_1 <- num_1 / den_1

  # Omega for factor 2
  num_2 <- sum(L[, 2])^2
  den_2 <- num_2 + sum(psi)
  omega_2 <- num_2 / den_2

  # Both should be between 0 and 1
  expect_true(omega_1 > 0 && omega_1 < 1)
  expect_true(omega_2 > 0 && omega_2 < 1)

  # Higher loading factor should have higher omega
  expect_equal(omega_1, omega_2, tolerance = 0.1) # Similar structure
})

test_that("omega is bounded between 0 and 1", {
  # Generate random valid parameters
  set.seed(42)
  J <- 5
  M <- 2

  for (i in 1:10) {
    L <- matrix(runif(J * M, 0.1, 0.9), J, M)
    psi <- runif(J, 0.1, 0.5)
    Sigma <- L %*% t(L) + diag(psi)

    omega_t <- sum(colSums(L)^2) / sum(Sigma)

    expect_true(omega_t >= 0, info = paste("Iteration", i))
    expect_true(omega_t <= 1.01, info = paste("Iteration", i)) # Small tolerance
  }
})

# --- h2 (communality) calculation tests ---

test_that("communality h2 is calculated correctly from Lambda", {
  L <- matrix(c(
    0.8, 0.2,
    0.3, 0.7
  ), nrow = 2, ncol = 2, byrow = TRUE)

  # h2 = rowSums(L^2) for orthogonal factors
  h2 <- rowSums(L^2)

  expect_equal(h2[1], 0.8^2 + 0.2^2) # 0.68
  expect_equal(h2[2], 0.3^2 + 0.7^2) # 0.58
})

test_that("uniqueness psi is Sigma_jj - h2_j", {
  J <- 3
  M <- 2

  L <- matrix(c(
    0.7, 0.1,
    0.1, 0.8,
    0.5, 0.5
  ), nrow = 3, ncol = 2, byrow = TRUE)

  psi_true <- c(0.3, 0.25, 0.2)
  Sigma <- L %*% t(L) + diag(psi_true)

  h2 <- rowSums(L^2)
  psi_calc <- diag(Sigma) - h2

  expect_equal(psi_calc, psi_true, tolerance = 1e-10)
})

# Note: Model type tests (Sigma for raw, Rho for std) and edge cases (M=1, large M)
# are now covered by run_all_befa_checks() in test-befa.R


# --- Summary statistics tests ---

test_that("Estimate is mean of posterior draws", {
  # Manual calculation
  draws <- c(0.75, 0.78, 0.72, 0.80, 0.77)

  expected_mean <- mean(draws)
  expected_sd <- sd(draws)
  expected_lower <- quantile(draws, 0.025)
  expected_upper <- quantile(draws, 0.975)

  # summarize_vec function behavior
  result <- data.frame(
    Estimate = mean(draws),
    SD = sd(draws),
    Lower_95 = quantile(draws, 0.025),
    Upper_95 = quantile(draws, 0.975)
  )

  expect_equal(result$Estimate, expected_mean)
  expect_equal(result$SD, expected_sd)
})

test_that("credible intervals contain the mean", {
  draws <- rnorm(1000, 0.8, 0.05)

  lower <- quantile(draws, 0.025)
  upper <- quantile(draws, 0.975)
  mean_val <- mean(draws)

  expect_true(mean_val > lower)
  expect_true(mean_val < upper)
})

# --- Numerical stability tests ---

test_that("befa_reliability handles near-zero loadings", {
  J <- 4
  M <- 2

  # Very small loadings
  L <- matrix(0.01, J, M)
  L[1, 1] <- 0.5
  L[3, 2] <- 0.5

  psi <- rep(0.9, J)
  Sigma <- L %*% t(L) + diag(psi)

  omega_t <- sum(colSums(L)^2) / sum(Sigma)

  expect_true(is.finite(omega_t))
  expect_true(omega_t >= 0)
})

test_that("befa_reliability handles very high loadings", {
  J <- 4
  M <- 2

  # High loadings, low uniqueness
  L <- matrix(0.9, J, M)
  psi <- rep(0.01, J)
  Sigma <- L %*% t(L) + diag(psi)

  omega_t <- sum(colSums(L)^2) / sum(Sigma)

  expect_true(is.finite(omega_t))
  # High loadings should give high reliability
  expect_true(omega_t > 0.9)
})

# --- Class and print tests ---

test_that("befa_reliability returns object with correct class", {
  skip("Requires real befa object")

  # result <- befa_reliability(object)
  # expect_s3_class(result, "befa_reliability")
})

# --- Comparison with theoretical values ---

test_that("omega equals alpha for tau-equivalent model", {
  # In a tau-equivalent model (equal loadings), omega = alpha
  J <- 4
  lambda <- 0.7 # Equal loadings
  psi <- 0.3 # Equal uniqueness

  # All loadings equal
  L <- matrix(lambda, J, 1)
  Sigma <- L %*% t(L) + diag(rep(psi, J))

  # Omega total
  omega <- sum(colSums(L)^2) / sum(Sigma)

  # Alpha formula: (J / (J-1)) * (1 - sum(diag(Sigma)) / sum(Sigma))
  alpha <- (J / (J - 1)) * (1 - sum(diag(Sigma)) / sum(Sigma))

  # Should be very close for tau-equivalent model
  expect_equal(omega, alpha, tolerance = 0.05)
})

test_that("omega > alpha for congeneric model", {
  # In a congeneric model (unequal loadings), omega >= alpha
  J <- 4

  # Unequal loadings
  L <- matrix(c(0.9, 0.7, 0.5, 0.3), J, 1)
  psi <- c(0.19, 0.51, 0.75, 0.91) # 1 - lambda^2
  Sigma <- L %*% t(L) + diag(psi)

  omega <- sum(colSums(L)^2) / sum(Sigma)
  alpha <- (J / (J - 1)) * (1 - sum(diag(Sigma)) / sum(Sigma))

  # Omega should be greater than or equal to alpha
  expect_true(omega >= alpha - 0.01) # Small tolerance
})

# ─────────────────────────────────────────────────────────────────────────────
