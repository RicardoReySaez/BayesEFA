# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: set_default_priors function
# ─────────────────────────────────────────────────────────────────────────────

test_that("get_default_priors returns lambda for uni-forced case (M=1, unit_vector)", {
  result <- get_default_priors("cor", "unit_vector", n_factors = 1)
  expect_true("lambda" %in% names(result))
  expect_equal(result$lambda, c(1, 1)) # Beta prior
  expect_false("xi" %in% names(result)) # No xi for uni
  expect_false("h2" %in% names(result)) # No h2 for uni
})

test_that("get_default_priors returns xi and h2 for unit_vector with M > 1", {
  result <- get_default_priors("cor", "unit_vector", n_factors = 3)
  expect_equal(result$xi, 100)
  expect_equal(result$h2, c(1, 1))
  expect_false("lambda" %in% names(result))
})

test_that("get_default_priors returns lambda and psi for normal prior", {
  result <- get_default_priors("raw", "normal", n_factors = 2)
  expect_equal(result$lambda, c(0, 10))
  expect_equal(result$psi, c(0.5, 0.5))
})

test_that("get_default_priors adds nu and sigma for raw model", {
  result <- get_default_priors("raw", "unit_vector", n_factors = 2)
  expect_equal(result$nu, c(0, 40))
  expect_equal(result$sigma, c(3, 0, 2.5))
})

test_that("get_default_priors does NOT add nu/sigma for std model", {
  result <- get_default_priors("cor", "unit_vector", n_factors = 2)
  expect_false("nu" %in% names(result))
  expect_false("sigma" %in% names(result))
})

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: prepare_befa_priors
# ─────────────────────────────────────────────────────────────────────────────

test_that("prepare_befa_priors returns defaults when prior_user is empty", {
  result <- prepare_befa_priors(list(), "cor", "unit_vector", n_factors = 2)
  expect_equal(result$xi, 100)
  expect_equal(result$h2, c(1, 1))
})

test_that("prepare_befa_priors merges user priors with defaults", {
  user_prior <- list(xi = 50)
  result <- prepare_befa_priors(user_prior, "cor", "unit_vector", n_factors = 2)
  expect_equal(result$xi, 50) # User override
  expect_equal(result$h2, c(1, 1)) # Default preserved
})

test_that("prepare_befa_priors warns on unknown prior parameters", {
  user_prior <- list(unknown_param = 999, another_bad = "foo")
  expect_warning(
    prepare_befa_priors(user_prior, "cor", "unit_vector", n_factors = 2),
    "not valid for this configuration"
  )
})

test_that("prepare_befa_priors errors if nu is wrong length (raw model)", {
  user_prior <- list(nu = c(0, 1, 2)) # Should be length 2
  expect_error(
    prepare_befa_priors(user_prior, "raw", "unit_vector", n_factors = 2),
    "vector of length 2"
  )
})

test_that("prepare_befa_priors errors if nu SD is non-positive", {
  user_prior <- list(nu = c(0, 0)) # SD = 0

  expect_error(
    prepare_befa_priors(user_prior, "raw", "unit_vector", n_factors = 2),
    "SD must be positive"
  )
})

test_that("prepare_befa_priors errors if sigma is wrong length", {
  user_prior <- list(sigma = c(3, 0)) # Should be length 3
  expect_error(
    prepare_befa_priors(user_prior, "raw", "unit_vector", n_factors = 2),
    "vector of length 3"
  )
})

test_that("prepare_befa_priors errors if lambda is wrong length (uni case)", {
  user_prior <- list(lambda = c(1)) # Should be length 2
  expect_error(
    prepare_befa_priors(user_prior, "cor", "unit_vector", n_factors = 1),
    "vector of length 2"
  )
})

test_that("prepare_befa_priors errors if lambda has non-positive Beta params (uni)", {
  user_prior <- list(lambda = c(0, 1)) # First param = 0
  expect_error(
    prepare_befa_priors(user_prior, "cor", "unit_vector", n_factors = 1),
    "must be positive"
  )
})

test_that("prepare_befa_priors errors if xi is not a single numeric", {
  user_prior <- list(xi = c(100, 200))
  expect_error(
    prepare_befa_priors(user_prior, "cor", "unit_vector", n_factors = 2),
    "single numeric"
  )
})

test_that("prepare_befa_priors errors if h2 is wrong length", {
  user_prior <- list(h2 = c(1)) # Should be length 2
  expect_error(
    prepare_befa_priors(user_prior, "cor", "unit_vector", n_factors = 2),
    "vector of length 2"
  )
})

test_that("prepare_befa_priors errors if lambda is wrong length (normal)", {
  user_prior <- list(lambda = c(0))
  expect_error(
    prepare_befa_priors(user_prior, "raw", "normal", n_factors = 2),
    "vector of length 2"
  )
})

test_that("prepare_befa_priors errors if psi is wrong length (normal)", {
  user_prior <- list(psi = c(0.5))
  expect_error(
    prepare_befa_priors(user_prior, "raw", "normal", n_factors = 2),
    "vector of length 2"
  )
})

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: append_priors_to_data
# ─────────────────────────────────────────────────────────────────────────────

# Helper: minimal stan_data structure
get_minimal_stan_data <- function(M = 2) {
  list(N = 100, J = 5, M = M, R_obs = diag(5))
}

test_that("append_priors_to_data adds pr_nu and pr_sigma for raw model", {
  stan_data <- get_minimal_stan_data()
  final_prior <- list(nu = c(0, 40), sigma = c(3, 0, 2.5), xi = 100, h2 = c(1, 1))
  result <- append_priors_to_data(stan_data, final_prior, "raw", "unit_vector")

  expect_equal(result$pr_nu, c(0, 40))
  expect_equal(result$pr_sigma, c(3, 0, 2.5))
})

test_that("append_priors_to_data does NOT add nu/sigma for std model", {
  stan_data <- get_minimal_stan_data()
  final_prior <- list(xi = 100, h2 = c(1, 1))
  result <- append_priors_to_data(stan_data, final_prior, "cor", "unit_vector")

  expect_equal(result$pr_nu, c(0, 10))
  expect_equal(result$pr_sigma, c(3, 0, 2.5))
})

test_that("append_priors_to_data adds pr_Lambda for uni-forced case (M=1)", {
  stan_data <- get_minimal_stan_data(M = 1)
  final_prior <- list(lambda = c(1, 1))
  result <- append_priors_to_data(stan_data, final_prior, "cor", "unit_vector")

  expect_equal(result$pr_Lambda, c(1, 1))
  expect_equal(result$pr_xi, 100)
})

test_that("append_priors_to_data adds pr_xi and pr_h2 for unit_vector (M > 1)", {
  stan_data <- get_minimal_stan_data(M = 3)
  final_prior <- list(xi = 100, h2 = c(1, 1))
  result <- append_priors_to_data(stan_data, final_prior, "cor", "unit_vector")

  expect_equal(result$pr_xi, 100)
  expect_equal(result$pr_h2, c(1, 1))
})

test_that("append_priors_to_data adds pr_Lambda and pr_psi for normal prior", {
  stan_data <- get_minimal_stan_data()
  final_prior <- list(lambda = c(0, 10), psi = c(0.5, 0.5), nu = c(0, 40), sigma = c(3, 0, 2.5))
  result <- append_priors_to_data(stan_data, final_prior, "raw", "normal")

  expect_equal(result$pr_Lambda, c(0, 10))
  expect_equal(result$pr_psi, c(0.5, 0.5))
})

test_that("append_priors_to_data preserves original stan_data fields", {
  stan_data <- get_minimal_stan_data()
  final_prior <- list(xi = 100, h2 = c(1, 1))
  result <- append_priors_to_data(stan_data, final_prior, "cor", "unit_vector")

  # Original fields preserved
  expect_equal(result$N, 100)
  expect_equal(result$J, 5)
  expect_equal(result$M, 2)
})

# ─────────────────────────────────────────────────────────────────────────────
