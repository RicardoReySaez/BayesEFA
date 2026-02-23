# ═══════════════════════════════════════════════════════════════════════════════
#                    Tests for FIML Utilities (utils_FIML.R)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("compute_fiml_loglik matches mvnfast with complete data", {
  set.seed(123)
  N <- 50
  J <- 4
  Mu <- c(1, 2, 3, 4)
  Sigma <- matrix(c(
    1.0, 0.5, 0.3, 0.2,
    0.5, 1.0, 0.4, 0.3,
    0.3, 0.4, 1.0, 0.5,
    0.2, 0.3, 0.5, 1.0
  ), nrow = J)

  Y <- mvnfast::rmvn(N, mu = Mu, sigma = Sigma)

  # Build pre-split data_fiml list (all complete, no missing)
  comp_idx <- seq_len(N)
  data_fiml <- list(
    comp_idx = comp_idx,
    miss_idx = integer(0),
    Y_comp   = Y,
    Y_miss   = Y[integer(0), , drop = FALSE]
  )

  # FIML result
  ll_fiml <- compute_fiml_loglik(data_fiml, mu = Mu, Sigma = Sigma)

  # Direct mvnfast result
  ll_direct <- mvnfast::dmvn(Y, mu = Mu, sigma = Sigma, log = TRUE)

  expect_equal(ll_fiml, ll_direct, tolerance = 1e-10)
})


test_that("compute_fiml_loglik handles NAs by marginalization", {
  set.seed(456)
  J <- 3
  Mu <- c(0, 0, 0)
  Sigma <- diag(1, J)

  # Single observation with one NA
  Y <- matrix(c(1.0, NA, 0.5), nrow = 1)

  # Build pre-split data_fiml list (all incomplete)
  data_fiml <- list(
    comp_idx = integer(0),
    miss_idx = 1L,
    Y_comp   = Y[integer(0), , drop = FALSE],
    Y_miss   = Y
  )

  ll_fiml <- compute_fiml_loglik(data_fiml, mu = Mu, Sigma = Sigma)

  # Manual: marginal N(0,1) for observed elements only
  ll_expected <- sum(dnorm(c(1.0, 0.5), mean = 0, sd = 1, log = TRUE))

  expect_equal(ll_fiml[1], ll_expected, tolerance = 1e-10)
})


test_that("compute_fiml_loglik returns 0 for fully missing row", {
  Y <- matrix(c(NA, NA, NA), nrow = 1)
  Mu <- c(0, 0, 0)
  Sigma <- diag(1, 3)

  # Build pre-split data_fiml list (all incomplete)
  data_fiml <- list(
    comp_idx = integer(0),
    miss_idx = 1L,
    Y_comp   = Y[integer(0), , drop = FALSE],
    Y_miss   = Y
  )

  ll <- compute_fiml_loglik(data_fiml, mu = Mu, Sigma = Sigma)
  expect_equal(ll[1], 0)
})


test_that("estimate_h1_moments_em recovers true moments with complete data", {
  set.seed(789)
  N <- 200
  J <- 3
  Mu_true <- c(2, 5, 8)
  Sigma_true <- matrix(c(
    1.0, 0.6, 0.3,
    0.6, 2.0, 0.5,
    0.3, 0.5, 1.5
  ), nrow = J)

  Y <- mvnfast::rmvn(N, mu = Mu_true, sigma = Sigma_true)

  # EM estimate
  em_est <- estimate_h1_moments_em(Y)

  # Should match sample moments exactly (no missing data)
  expect_equal(em_est$Mu, colMeans(Y), tolerance = 1e-6)
  expect_equal(em_est$Sigma, cov(Y) * (N - 1) / N, tolerance = 1e-6)
})


test_that("estimate_h1_moments_em converges with missing data", {
  set.seed(101112)
  N <- 100
  J <- 3
  Mu_true <- c(0, 0, 0)
  Sigma_true <- diag(1, J)

  Y <- mvnfast::rmvn(N, mu = Mu_true, sigma = Sigma_true)

  # Introduce ~20% missing at random
  na_mask <- matrix(runif(N * J) < 0.2, nrow = N, ncol = J)
  Y[na_mask] <- NA

  # EM should converge without error
  em_est <- estimate_h1_moments_em(Y)

  # Estimates should be reasonable (close to true values)
  expect_true(all(abs(em_est$Mu) < 0.5)) # Near 0
  expect_true(all(diag(em_est$Sigma) > 0.5 & diag(em_est$Sigma) < 1.5)) # Near 1
})


test_that("em_estep returns correct dimensions", {
  N <- 20
  J <- 4
  Y <- matrix(rnorm(N * J), nrow = N)
  Mu <- rep(0, J)
  Sigma <- diag(1, J)

  estep <- em_E_step(Y, Mu, Sigma)

  expect_length(estep$T1, J)
  expect_equal(dim(estep$T2), c(J, J))
})
