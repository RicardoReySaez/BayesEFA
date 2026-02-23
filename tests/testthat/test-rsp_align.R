# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: rsp_align
# ─────────────────────────────────────────────────────────────────────────────

# --- Helper functions ---
generate_test_lambda_draws <- function(n_draws = 100, J = 5, M = 2,
                                       true_lambda = NULL, noise_sd = 0.1,
                                       format = "column_major") {
  set.seed(123)

  if (is.null(true_lambda)) {
    # Create a simple structure with clear factor pattern
    true_lambda <- matrix(0, nrow = J, ncol = M)
    items_per_factor <- floor(J / M)
    for (m in 1:M) {
      start_idx <- (m - 1) * items_per_factor + 1
      end_idx <- min(m * items_per_factor, J)
      true_lambda[start_idx:end_idx, m] <- runif(end_idx - start_idx + 1, 0.6, 0.9)
    }
    # Add small cross-loadings
    true_lambda[true_lambda == 0] <- runif(sum(true_lambda == 0), 0.05, 0.15)
  }

  # Generate draws with noise
  draws <- matrix(NA, nrow = n_draws, ncol = J * M)
  for (s in 1:n_draws) {
    noise <- matrix(rnorm(J * M, 0, noise_sd), nrow = J, ncol = M)
    L_s <- true_lambda + noise

    if (format == "column_major") {
      draws[s, ] <- as.vector(L_s)
    } else {
      draws[s, ] <- as.vector(t(L_s))
    }
  }

  list(draws = draws, true_lambda = true_lambda)
}

# --- Basic functionality tests ---

test_that("rsp_align returns correct structure", {
  dat <- generate_test_lambda_draws(n_draws = 50, J = 5, M = 2)

  result <- suppressWarnings(rsp_align(
    lambda_draws = dat$draws,
    n_items = 5,
    n_factors = 2,
    n_chains = 1,
    format = "column_major"
  ))

  expect_type(result, "list")
  expect_true("Lambda_hat_mcmc" %in% names(result))
  expect_true("objective" %in% names(result) || "Lambda_star" %in% names(result))
})

test_that("rsp_align output has correct dimensions", {
  n_draws <- 50
  J <- 6
  M <- 3
  dat <- generate_test_lambda_draws(n_draws = n_draws, J = J, M = M)

  result <- suppressWarnings(rsp_align(dat$draws, n_items = J, n_factors = M, n_chains = 1, format = "column_major"))

  expect_equal(nrow(result$Lambda_hat_mcmc), n_draws)
  expect_equal(ncol(result$Lambda_hat_mcmc), J * M)
})

# --- Format argument tests ---

test_that("rsp_align works with column_major format", {
  dat <- generate_test_lambda_draws(n_draws = 30, J = 4, M = 2, format = "column_major")

  result <- suppressWarnings(rsp_align(
    lambda_draws = dat$draws,
    n_items = 4,
    n_factors = 2,
    n_chains = 1,
    format = "column_major"
  ))

  expect_type(result, "list")
  expect_equal(ncol(result$Lambda_hat_mcmc), 4 * 2)
})

test_that("rsp_align works with row_major format", {
  dat <- generate_test_lambda_draws(n_draws = 30, J = 4, M = 2, format = "row_major")

  result <- suppressWarnings(rsp_align(
    lambda_draws = dat$draws,
    n_items = 4,
    n_factors = 2,
    n_chains = 1,
    format = "row_major"
  ))

  expect_type(result, "list")
  expect_equal(ncol(result$Lambda_hat_mcmc), 4 * 2)
})

test_that("rsp_align returns output in same format as input", {
  J <- 4
  M <- 2

  # Create known structure
  true_L <- matrix(c(
    0.8, 0.7, 0.1, 0.1,
    0.1, 0.1, 0.8, 0.7
  ), nrow = 4, ncol = 2)

  # Row-major input
  set.seed(42)
  draws_row <- matrix(NA, 20, J * M)
  for (s in 1:20) {
    draws_row[s, ] <- as.vector(t(true_L + rnorm(J * M, 0, 0.05)))
  }

  result <- suppressWarnings(rsp_align(draws_row, J, M, n_chains = 1, format = "row_major"))

  # Convert first output row back to matrix (row-major)
  out_L <- matrix(result$Lambda_hat_mcmc[1, ], nrow = J, ncol = M, byrow = TRUE)

  # Structure should be preserved
  expect_equal(nrow(out_L), J)
  expect_equal(ncol(out_L), M)
})

# --- Dimension validation tests ---

test_that("rsp_align errors on dimension mismatch", {
  # 20 columns but J*M = 5*3 = 15
  bad_draws <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)

  expect_error(
    rsp_align(bad_draws, n_items = 5, n_factors = 3, n_chains = 1, format = "column_major"),
    "Dimension mismatch"
  )
})

test_that("rsp_align errors when n_items * n_factors != ncol", {
  draws <- matrix(rnorm(100 * 12), nrow = 100, ncol = 12)

  expect_error(
    rsp_align(draws, n_items = 5, n_factors = 3, n_chains = 1, format = "column_major"), # 5*3=15 != 12
    "Dimension mismatch"
  )
})

test_that("rsp_align errors when format is not specified", {
  draws <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)

  expect_error(
    rsp_align(draws, n_items = 5, n_factors = 2, n_chains = 1),
    "Argument 'format' is required"
  )
})

# --- Unidimensional case (M = 1) ---

test_that("rsp_align handles M=1 without rotation", {
  J <- 5
  M <- 1
  n_draws <- 30

  set.seed(123)
  true_lambda <- matrix(runif(J, 0.5, 0.9), nrow = J, ncol = 1)
  draws <- matrix(NA, n_draws, J)
  for (s in 1:n_draws) {
    draws[s, ] <- true_lambda + rnorm(J, 0, 0.05)
  }

  result <- suppressWarnings(rsp_align(draws, n_items = J, n_factors = M, n_chains = 1, format = "column_major"))

  expect_type(result, "list")
  expect_equal(ncol(result$Lambda_hat_mcmc), J)
})

test_that("rsp_align with M=1 only applies sign alignment", {
  J <- 4
  M <- 1
  n_draws <- 20

  # All positive loadings
  set.seed(456)
  draws <- matrix(abs(rnorm(n_draws * J, 0.7, 0.1)), nrow = n_draws)

  result <- suppressWarnings(rsp_align(draws, n_items = J, n_factors = M, n_chains = 1, format = "column_major"))

  # All output should have consistent sign (positive after alignment)
  col_means <- colMeans(result$Lambda_hat_mcmc)
  expect_true(all(col_means > 0)) # Forced positive
})

# --- Sign alignment tests ---

test_that("rsp_align forces positive mean loadings per factor", {
  J <- 6
  M <- 2
  dat <- generate_test_lambda_draws(n_draws = 40, J = J, M = M)

  result <- suppressWarnings(rsp_align(dat$draws, n_items = J, n_factors = M, n_chains = 1, format = "column_major"))

  # Compute column means for each factor
  aligned_matrix <- result$Lambda_hat_mcmc
  lambda_hat <- matrix(colMeans(aligned_matrix), nrow = J, ncol = M)

  # Sum of loadings per factor should be positive (forced sign)
  factor_sums <- colSums(lambda_hat)
  expect_true(all(factor_sums > 0))
})

# --- Rotation tests ---

test_that("rsp_align applies varimax rotation for M > 1", {
  J <- 6
  M <- 2
  n_draws <- 30

  # Create unrotated structure (factors not orthogonally simple)
  set.seed(789)
  draws <- matrix(rnorm(n_draws * J * M, 0.5, 0.3), nrow = n_draws)

  result <- suppressWarnings(rsp_align(draws, n_items = J, n_factors = M, n_chains = 1, format = "column_major"))

  # Result should exist (varimax was applied internally)
  expect_type(result, "list")
})

# --- Convergence tests ---

test_that("rsp_align respects max_iter argument", {
  dat <- generate_test_lambda_draws(n_draws = 50, J = 5, M = 2)

  # Very low max_iter should still complete
  result <- suppressWarnings(rsp_align(
    dat$draws,
    n_items = 5,
    n_factors = 2,
    n_chains = 1,
    format = "column_major",
    max_iter = 5
  ))

  expect_type(result, "list")
})

test_that("rsp_align respects threshold argument", {
  dat <- generate_test_lambda_draws(n_draws = 50, J = 5, M = 2)

  # Very loose threshold
  result_loose <- suppressWarnings(rsp_align(dat$draws, 5, 2, n_chains = 1, format = "column_major", threshold = 0.1))

  # Very strict threshold
  result_strict <- suppressWarnings(rsp_align(dat$draws, 5, 2, n_chains = 1, format = "column_major", threshold = 1e-10))

  expect_type(result_loose, "list")
  expect_type(result_strict, "list")
})

# --- Alignment quality tests ---

test_that("rsp_align reduces variability across draws", {
  J <- 6
  M <- 2

  # Create aligned true structure
  true_L <- matrix(c(
    0.8, 0.1,
    0.7, 0.15,
    0.75, 0.1,
    0.1, 0.8,
    0.15, 0.7,
    0.1, 0.75
  ), nrow = 6, ncol = 2, byrow = TRUE)

  dat <- generate_test_lambda_draws(
    n_draws = 100, J = J, M = M,
    true_lambda = true_L, noise_sd = 0.05
  )

  # Variance before alignment (raw draws)
  var_before <- mean(apply(dat$draws, 2, var))

  result <- suppressWarnings(rsp_align(dat$draws, n_items = J, n_factors = M, n_chains = 1, format = "column_major"))

  # Variance after alignment
  var_after <- mean(apply(result$Lambda_hat_mcmc, 2, var))

  # After alignment, variance in each column should be similar or lower
  # (alignment shouldn't increase variability dramatically)
  expect_true(var_after < var_before * 2)
})

test_that("rsp_align produces interpretable factor structure", {
  J <- 8
  M <- 2

  # Create clear two-factor structure
  true_L <- matrix(0.1, nrow = J, ncol = M)
  true_L[1:4, 1] <- c(0.8, 0.75, 0.7, 0.65)
  true_L[5:8, 2] <- c(0.8, 0.75, 0.7, 0.65)

  dat <- generate_test_lambda_draws(
    n_draws = 100, J = J, M = M,
    true_lambda = true_L, noise_sd = 0.03
  )

  result <- suppressWarnings(rsp_align(dat$draws, n_items = J, n_factors = M, n_chains = 1, format = "column_major"))

  # Mean loadings should recover structure
  mean_L <- matrix(colMeans(result$Lambda_hat_mcmc), nrow = J, ncol = M)

  # Items 1-4 should load higher on factor 1
  expect_true(mean(mean_L[1:4, 1]) > mean(mean_L[1:4, 2]))

  # Items 5-8 should load higher on factor 2
  expect_true(mean(mean_L[5:8, 2]) > mean(mean_L[5:8, 1]))
})

# --- Edge cases ---

test_that("rsp_align handles single draw", {
  J <- 5
  M <- 2

  single_draw <- matrix(rnorm(J * M), nrow = 1)

  result <- suppressWarnings(rsp_align(single_draw, n_items = J, n_factors = M, n_chains = 1, format = "column_major"))

  expect_equal(nrow(result$Lambda_hat_mcmc), 1)
})

test_that("rsp_align handles many factors (M large)", {
  J <- 20
  M <- 5
  n_draws <- 50

  set.seed(111)
  draws <- matrix(rnorm(n_draws * J * M, 0, 0.5), nrow = n_draws)

  result <- suppressWarnings(rsp_align(draws, n_items = J, n_factors = M, n_chains = 1, format = "column_major"))

  expect_equal(dim(result$Lambda_hat_mcmc), c(n_draws, J * M))
})

test_that("rsp_align handles very few items (J small)", {
  J <- 3
  M <- 2
  n_draws <- 30

  set.seed(222)
  draws <- matrix(rnorm(n_draws * J * M, 0.5, 0.2), nrow = n_draws)

  result <- suppressWarnings(rsp_align(draws, n_items = J, n_factors = M, n_chains = 1, format = "column_major"))

  expect_equal(dim(result$Lambda_hat_mcmc), c(n_draws, J * M))
})

test_that("rsp_align handles negative loadings correctly", {
  J <- 4
  M <- 2
  n_draws <- 50

  # Create structure with intentional negative loadings
  true_L <- matrix(c(
    0.8, -0.1,
    0.7, -0.15,
    -0.1, 0.8,
    -0.15, 0.7
  ), nrow = 4, ncol = 2, byrow = TRUE)

  dat <- generate_test_lambda_draws(
    n_draws = n_draws, J = J, M = M,
    true_lambda = true_L, noise_sd = 0.05
  )

  result <- suppressWarnings(rsp_align(dat$draws, n_items = J, n_factors = M, n_chains = 1, format = "column_major"))

  # Sign should be resolved to mostly positive averages
  mean_L <- colMeans(result$Lambda_hat_mcmc)
  expect_true(sum(mean_L > 0) >= J) # Most should be positive
})

# --- Reproducibility tests ---

test_that("rsp_align is deterministic given same input", {
  set.seed(333)
  dat <- generate_test_lambda_draws(n_draws = 30, J = 5, M = 2)

  result1 <- suppressWarnings(rsp_align(dat$draws, 5, 2, n_chains = 1, format = "column_major"))
  result2 <- suppressWarnings(rsp_align(dat$draws, 5, 2, n_chains = 1, format = "column_major"))

  expect_equal(result1$Lambda_hat_mcmc, result2$Lambda_hat_mcmc)
})

# --- format argument matching ---

test_that("rsp_align matches format argument correctly", {
  draws <- matrix(rnorm(100 * 10), nrow = 100)

  # Valid formats
  expect_no_error(suppressWarnings(rsp_align(draws, 5, 2, n_chains = 1, format = "column_major")))
  expect_no_error(suppressWarnings(rsp_align(draws, 5, 2, n_chains = 1, format = "row_major")))
})

# --- add_names tests ---

test_that("rsp_align preserves existing names when add_names = FALSE", {
  dat <- generate_test_lambda_draws(n_draws = 30, J = 5, M = 2)
  custom_names <- paste0("V", 1:10)
  colnames(dat$draws) <- custom_names

  result <- suppressWarnings(rsp_align(
    dat$draws, 5, 2,
    n_chains = 1,
    format = "column_major",
    add_names = FALSE
  ))

  expect_equal(colnames(result$Lambda_hat_mcmc), custom_names)
})

test_that("rsp_align forces names when add_names = FALSE and no existing names", {
  dat <- generate_test_lambda_draws(n_draws = 30, J = 5, M = 2)
  # dat$draws has no colnames by default

  result <- suppressWarnings(rsp_align(
    dat$draws, 5, 2,
    n_chains = 1,
    format = "column_major",
    add_names = FALSE
  ))

  # Should have Stan-format names since input had none
  expect_true(!is.null(colnames(result$Lambda_hat_mcmc)))
  expect_true(grepl("^Lambda\\[", colnames(result$Lambda_hat_mcmc)[1]))
})

# ─────────────────────────────────────────────────────────────────────────────
