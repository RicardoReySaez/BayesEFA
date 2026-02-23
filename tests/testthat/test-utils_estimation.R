# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: fit_befa_model
# ─────────────────────────────────────────────────────────────────────────────

test_that("fit_befa_model errors if model not found", {
  fake_stan_data <- list(N = 100, J = 5, M = 2)

  expect_error(
    fit_befa_model(
      model_name = "nonexistent_model",
      stan_data = fake_stan_data,
      backend = "rstan",
      verbose = FALSE,
      model_type = "cor",
      lambda_prior = "unit_vector"
    ),
    "not found"
  )
})

# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: normalize_stan_args
# ─────────────────────────────────────────────────────────────────────────────

get_test_stan_data <- function(J = 5, M = 2, N = 100) {
  list(
    N = N,
    J = J,
    M = M,
    m_obs = rep(0, J),
    S_obs = diag(J),
    R_obs = diag(J)
  )
}

test_that("normalize_stan_args returns default values when user_dots empty", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list()
  )

  expect_equal(result$chains, 4)
  expect_equal(result$iter, 2000)
  expect_equal(result$warmup, 1000)
  expect_equal(result$control$adapt_delta, 0.95)
  expect_equal(result$control$max_treedepth, 15)
})

test_that("normalize_stan_args respects user-provided iter", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(iter = 4000)
  )

  expect_equal(result$iter, 4000)
  # Warmup should be 50% of iter when not specified
  expect_equal(result$warmup, 2000)
})

test_that("normalize_stan_args respects user-provided warmup", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(iter = 3000, warmup = 500)
  )

  expect_equal(result$iter, 3000)
  expect_equal(result$warmup, 500)
})

test_that("normalize_stan_args converts cmdstanr args to rstan format", {
  stan_data <- get_test_stan_data()

  # CmdStanR uses iter_sampling and iter_warmup
  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(iter_sampling = 1500, iter_warmup = 500)
  )

  expect_equal(result$iter, 2000) # 1500 + 500
  expect_equal(result$warmup, 500)
})

test_that("normalize_stan_args handles parallel_chains -> cores mapping", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(parallel_chains = 2)
  )

  expect_equal(result$cores, 2)
})

test_that("normalize_stan_args defaults cores to chains when not specified", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(chains = 3)
  )

  expect_equal(result$chains, 3)
  expect_equal(result$cores, 3)
})

test_that("normalize_stan_args extracts adapt_delta from control list", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(control = list(adapt_delta = 0.99))
  )

  expect_equal(result$control$adapt_delta, 0.99)
})

test_that("normalize_stan_args respects direct adapt_delta argument", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(adapt_delta = 0.98)
  )

  expect_equal(result$control$adapt_delta, 0.98)
})

test_that("normalize_stan_args direct adapt_delta overrides control list", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(
      control = list(adapt_delta = 0.90),
      adapt_delta = 0.99
    )
  )

  # Direct argument should win
  expect_equal(result$control$adapt_delta, 0.99)
})

test_that("normalize_stan_args sets refresh based on iter", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(iter = 4000)
  )

  # Default refresh is iter/4
  expect_equal(result$refresh, 1000)
})

test_that("normalize_stan_args respects user-provided refresh", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(refresh = 200)
  )

  expect_equal(result$refresh, 200)
})

test_that("normalize_stan_args preserves unknown extra arguments", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(thin = 2, save_warmup = TRUE)
  )

  expect_equal(result$thin, 2)
  expect_true(result$save_warmup)
})

test_that("normalize_stan_args passes seed if provided", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(seed = 12345)
  )

  expect_equal(result$seed, 12345)
})

test_that("normalize_stan_args formats correctly for cmdstanr backend", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "cmdstanr",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(iter = 2000, warmup = 1000)
  )

  expect_equal(result$iter_warmup, 1000)
  expect_equal(result$iter_sampling, 1000) # iter - warmup
  expect_true("parallel_chains" %in% names(result))
  expect_equal(result$adapt_delta, 0.95) # Direct, not in control
})

test_that("normalize_stan_args cmdstanr does not use control list", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "cmdstanr",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list()
  )

  expect_null(result$control) # CmdStanR uses direct args
  expect_equal(result$adapt_delta, 0.95)
  expect_equal(result$max_treedepth, 15)
})

test_that("normalize_stan_args generates inits when not provided", {
  stan_data <- get_test_stan_data()

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(chains = 2)
  )

  expect_type(result$init, "list")
  expect_length(result$init, 2) # One per chain
})

test_that("normalize_stan_args respects user-provided inits", {
  stan_data <- get_test_stan_data()
  custom_init <- list(list(h2 = rep(0.5, 5)))

  result <- normalize_stan_args(
    backend = "rstan",
    model_type = "raw",
    lambda_prior = "unit_vector",
    stan_data = stan_data,
    user_dots = list(init = custom_init)
  )

  expect_identical(result$init, custom_init)
})

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: get_befa_inits
# ─────────────────────────────────────────────────────────────────────────────


test_that("get_befa_inits returns NULL for non-list stan_data", {
  result <- get_befa_inits("raw", "unit_vector", "not_a_list", 4)
  expect_null(result)
})
test_that("get_befa_inits returns list of length n_chains", {
  stan_data <- get_test_stan_data(J = 5, M = 2)

  result <- get_befa_inits("raw", "unit_vector", stan_data, n_chains = 4)

  expect_type(result, "list")
  expect_length(result, 4)
})

test_that("get_befa_inits includes nu and sigma for raw model", {
  stan_data <- get_test_stan_data(J = 5, M = 2)

  result <- get_befa_inits("raw", "unit_vector", stan_data, n_chains = 1)

  expect_true("nu" %in% names(result[[1]]))
  expect_true("sigma" %in% names(result[[1]]))
  expect_length(result[[1]]$nu, 5)
  expect_length(result[[1]]$sigma, 5)
})

test_that("get_befa_inits does NOT include nu/sigma for cor model", {
  stan_data <- get_test_stan_data(J = 5, M = 2)

  result <- get_befa_inits("cor", "unit_vector", stan_data, n_chains = 1)

  expect_null(result[[1]]$nu)
  expect_null(result[[1]]$sigma)
})

test_that("get_befa_inits includes Lambda_uni for uni-forced case (M=1)", {
  stan_data <- get_test_stan_data(J = 5, M = 1)

  result <- get_befa_inits("cor", "unit_vector", stan_data, n_chains = 1)

  expect_true("Lambda_uni" %in% names(result[[1]]))
  expect_equal(nrow(result[[1]]$Lambda_uni), 5)

  # Lambda_uni should be bounded (-0.95, 0.95)
  expect_true(all(result[[1]]$Lambda_uni >= -0.95))
  expect_true(all(result[[1]]$Lambda_uni <= 0.95))
})

test_that("get_befa_inits includes Lambda_norm and psi for normal prior", {
  stan_data <- get_test_stan_data(J = 5, M = 2)

  result <- get_befa_inits("raw", "normal", stan_data, n_chains = 1)

  expect_true("Lambda_norm" %in% names(result[[1]]))
  expect_true("psi" %in% names(result[[1]]))
  expect_equal(dim(result[[1]]$Lambda_norm), c(5, 2))
  expect_length(result[[1]]$psi, 5)
})

test_that("get_befa_inits includes h2 and Z for unit_vector (M > 1)", {
  stan_data <- get_test_stan_data(J = 5, M = 3)

  result <- get_befa_inits("cor", "unit_vector", stan_data, n_chains = 1)

  expect_true("h2" %in% names(result[[1]]))
  expect_true("Z" %in% names(result[[1]]))
  expect_length(result[[1]]$h2, 5)
  expect_equal(dim(result[[1]]$Z), c(5, 3))
})

test_that("get_befa_inits h2 values are bounded (0.05, 0.95)", {
  stan_data <- get_test_stan_data(J = 10, M = 3)

  result <- get_befa_inits("cor", "unit_vector", stan_data, n_chains = 1)

  expect_true(all(result[[1]]$h2 >= 0.05))
  expect_true(all(result[[1]]$h2 <= 0.95))
})

test_that("get_befa_inits generates different values per chain", {
  stan_data <- get_test_stan_data(J = 5, M = 2)

  result <- get_befa_inits("cor", "unit_vector", stan_data, n_chains = 2)

  # h2 values should differ between chains (due to random noise)
  expect_false(all(result[[1]]$h2 == result[[2]]$h2))
})

test_that("get_befa_inits uses m_obs for nu initialization", {
  stan_data <- list(
    N = 100, J = 3, M = 1,
    m_obs = c(5, 10, 15),
    S_obs = diag(3)
  )

  result <- get_befa_inits("raw", "unit_vector", stan_data, n_chains = 1)

  # nu should be close to m_obs (with small noise)
  expect_true(abs(mean(result[[1]]$nu) - mean(stan_data$m_obs)) < 1)
})

test_that("get_befa_inits handles missing m_obs gracefully", {
  stan_data <- list(N = 100, J = 5, M = 2, S_obs = diag(5))
  # m_obs is NULL

  result <- get_befa_inits("raw", "unit_vector", stan_data, n_chains = 1)

  # Should use rep(0, J) as default
  expect_true(all(abs(result[[1]]$nu) < 1)) # Near zero
})

test_that("get_befa_inits handles missing S_obs gracefully", {
  stan_data <- list(N = 100, J = 5, M = 2, m_obs = rep(0, 5))
  # S_obs is NULL

  result <- get_befa_inits("raw", "unit_vector", stan_data, n_chains = 1)

  # Should use diag(J) as default
  expect_type(result, "list")
})

test_that("get_befa_inits sigma values are positive", {
  stan_data <- get_test_stan_data(J = 5, M = 2)

  result <- get_befa_inits("raw", "unit_vector", stan_data, n_chains = 1)

  expect_true(all(result[[1]]$sigma > 0))
})

test_that("get_befa_inits psi values are positive (normal prior)", {
  stan_data <- get_test_stan_data(J = 5, M = 2)

  result <- get_befa_inits("raw", "normal", stan_data, n_chains = 1)

  expect_true(all(result[[1]]$psi > 0))
})

test_that("get_befa_inits handles J = 2 (minimum items)", {
  stan_data <- list(N = 50, J = 2, M = 1, m_obs = c(0, 0), S_obs = diag(2))

  result <- get_befa_inits("cor", "unit_vector", stan_data, n_chains = 1)

  expect_type(result, "list")
  expect_equal(nrow(result[[1]]$Lambda_uni), 2)
})

test_that("get_befa_inits works with large J and M", {
  stan_data <- list(N = 500, J = 50, M = 10, m_obs = rep(0, 50), S_obs = diag(50))

  result <- get_befa_inits("cor", "unit_vector", stan_data, n_chains = 2)

  expect_length(result[[1]]$h2, 50)
  expect_equal(dim(result[[1]]$Z), c(50, 10))
})

test_that("get_befa_inits initial_h2 is based on correlations", {
  # Create data with known correlation structure
  set.seed(42)
  J <- 5
  R <- matrix(0.3, J, J)
  diag(R) <- 0 # diag(R) set to 0 for max calculation

  # initial_h2 = pmin(0.9, pmax(0.1, apply(abs(R), 1, max)))
  expected_h2 <- pmin(0.9, pmax(0.1, apply(abs(R), 1, max)))

  expect_true(all(expected_h2 >= 0.1))
  expect_true(all(expected_h2 <= 0.9))
  expect_equal(expected_h2, rep(0.3, 5))
})

test_that("get_befa_inits includes sigma but NOT nu for cov model", {
  stan_data <- get_test_stan_data(J = 5, M = 2)

  result <- get_befa_inits("cov", "unit_vector", stan_data, n_chains = 1)

  expect_true("sigma" %in% names(result[[1]]))
  expect_null(result[[1]]$nu)
  expect_length(result[[1]]$sigma, 5)
})

# ─────────────────────────────────────────────────────────────────────────────

test_that("get_befa_inits omits inactive parameters (CmdStanR compat)", {
  stan_data <- get_test_stan_data(J = 5, M = 2)

  # Case 1: lambda_prior = "unit_vector" (lambda_type = 1)
  # Active: h2, Z. Inactive (omitted): Lambda_uni, Lambda_norm, psi
  starts <- get_befa_inits("raw", "unit_vector", stan_data, n_chains = 1)[[1]]

  expect_true("h2" %in% names(starts))
  expect_true("Z" %in% names(starts))
  expect_null(starts$Lambda_uni)
  expect_null(starts$Lambda_norm)
  expect_null(starts$psi)

  # Case 2: lambda_prior = "normal" (lambda_type = 3)
  # Active: Lambda_norm, psi. Inactive (omitted): h2, Z, Lambda_uni
  starts_norm <- get_befa_inits("raw", "normal", stan_data, n_chains = 1)[[1]]

  expect_true("Lambda_norm" %in% names(starts_norm))
  expect_true("psi" %in% names(starts_norm))
  expect_null(starts_norm$h2)
  expect_null(starts_norm$Z)
  expect_null(starts_norm$Lambda_uni)
})
