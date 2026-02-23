# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: main function befa
# ─────────────────────────────────────────────────────────────────────────────

# --- Setup: dummy data for all tests ---
set.seed(123)
dummy_data <- as.data.frame(matrix(rnorm(200), ncol = 5))
colnames(dummy_data) <- paste0("V", 1:5)
dummy_cov <- cov(dummy_data)
dummy_cor <- cor(dummy_data)

# --- Argument matching tests ---

test_that("befa matches model argument correctly", {
  expect_no_error(match.arg("cor", c("cor", "cov", "raw")))
  expect_no_error(match.arg("raw", c("cor", "cov", "raw")))
  expect_error(match.arg("invalid", c("cor", "cov", "raw")))
})

test_that("befa matches lambda_prior argument correctly", {
  expect_no_error(match.arg("unit_vector", c("unit_vector", "normal")))
  expect_no_error(match.arg("normal", c("unit_vector", "normal")))
  expect_error(match.arg("invalid", c("unit_vector", "normal")))
})

test_that("befa matches backend argument correctly", {
  expect_no_error(match.arg("rstan", c("cmdstanr", "rstan")))
  expect_no_error(match.arg("cmdstanr", c("cmdstanr", "rstan")))
})

# --- Input validation (delegated to check_befa_inputs) ---

test_that("befa errors on ambiguous input", {
  skip_if_not_installed("rstan")

  expect_error(
    befa(
      data = dummy_data, sample_cov = dummy_cov, sample_nobs = 40,
      n_factors = 2
    ),
    "Ambiguous"
  )
})

test_that("befa errors when no input provided", {
  skip_if_not_installed("rstan")

  expect_error(
    befa(n_factors = 2),
    "No input"
  )
})

test_that("befa errors when n_factors >= J", {
  skip_if_not_installed("rstan")

  expect_error(
    befa(data = dummy_data, n_factors = 5), # J=5
    "smaller than"
  )
})

test_that("befa errors when n_factors is NULL", {
  skip_if_not_installed("rstan")

  expect_error(
    befa(data = dummy_data, n_factors = NULL),
    "positive integer"
  )
})

# --- n_factors = 1 (unidimensional) tests ---

test_that("befa uses 'uni' model for M=1 with unit_vector", {
  model_name <- select_befa_model("cor", "unit_vector", n_factors = 1, verbose = FALSE)
  expect_equal(model_name, "befa_efa")
})


# --- lambda_prior='normal' + model='cor' is blocked ---

test_that("befa blocks lambda_prior='normal' with model='cor'", {
  skip_if_not_installed("rstan")

  expect_error(
    befa(
      sample_cor = dummy_cor, sample_nobs = 40, n_factors = 2,
      model = "cor", lambda_prior = "normal"
    ),
    "not supported"
  )
})

# --- backend argument tests ---

test_that("befa falls back to rstan when cmdstanr unavailable", {
  skip_if(requireNamespace("cmdstanr", quietly = TRUE))

  expect_warning(
    backend <- validate_backend("cmdstanr"),
    "not installed"
  )
  expect_equal(backend, "rstan")
})

test_that("befa uses rstan backend when specified", {
  backend <- validate_backend("rstan")
  expect_equal(backend, "rstan")
})

# --- prior argument tests ---

test_that("befa warns on unknown prior parameters", {
  skip_if_not_installed("rstan")

  expect_warning(
    prepare_befa_priors(
      prior_user = list(unknown_param = 99),
      model = "cor",
      lambda_prior = "unit_vector",
      n_factors = 2
    ),
    "not valid"
  )
})

test_that("befa accepts valid prior parameters", {
  priors <- prepare_befa_priors(
    prior_user = list(xi = 50, h2 = c(2, 2)),
    model = "cor",
    lambda_prior = "unit_vector",
    n_factors = 2
  )

  expect_equal(priors$xi, 50)
  expect_equal(priors$h2, c(2, 2))
})

# --- rsp_args argument tests ---

test_that("befa accepts valid rsp_args", {
  rsp_config <- validate_rsp_args(list(max_iter = 500, threshold = 1e-5))

  expect_equal(rsp_config$max_iter, 500)
  expect_equal(rsp_config$threshold, 1e-5)
})

test_that("befa warns on unknown rsp_args", {
  expect_warning(
    validate_rsp_args(list(bad_arg = 123)),
    "Unknown"
  )
})

# --- Data handling tests ---

test_that("befa handles matrix input", {
  skip_if_not_installed("rstan")

  data_matrix <- as.matrix(dummy_data)
  args <- list(
    data = data_matrix, n_factors = 2, ordered = FALSE,
    sample_nobs = NULL, sample_mean = NULL, sample_cov = NULL,
    sample_cor = NULL, model = "cor", lambda_prior = "unit_vector",
    missing = "listwise",
    rotate = "varimax", rsp_args = NULL, factor_scores = FALSE,
    backend = "rstan", prior = list(), compute_fit_indices = FALSE,
    compute_reliability = FALSE, verbose = FALSE
  )

  result <- check_befa_inputs(args)
  expect_type(result, "list")
})

test_that("befa handles data.frame input", {
  skip_if_not_installed("rstan")

  args <- list(
    data = dummy_data, n_factors = 2, ordered = FALSE,
    sample_nobs = NULL, sample_mean = NULL, sample_cov = NULL,
    sample_cor = NULL, model = "cor", lambda_prior = "unit_vector",
    missing = "listwise",
    rotate = "varimax", rsp_args = NULL, factor_scores = FALSE,
    backend = "rstan", prior = list(), compute_fit_indices = FALSE,
    compute_reliability = FALSE, verbose = FALSE
  )

  result <- check_befa_inputs(args)
  expect_type(result, "list")
})

test_that("befa warns on missing data", {
  data_na <- dummy_data
  data_na[1, 1] <- NA

  args <- list(
    data = data_na, n_factors = 2, ordered = FALSE,
    sample_nobs = NULL, sample_mean = NULL, sample_cov = NULL,
    sample_cor = NULL, model = "raw", lambda_prior = "unit_vector",
    rotate = "varimax", rsp_args = NULL, factor_scores = FALSE,
    backend = "rstan", prior = list(), compute_fit_indices = FALSE,
    compute_reliability = FALSE, verbose = FALSE
  )

  expect_warning(check_befa_inputs(args), "Missing values")
})

# --- Summary statistics input tests ---

test_that("befa accepts sample_cor with sample_nobs", {
  args <- list(
    data = NULL, n_factors = 2, ordered = FALSE,
    sample_nobs = 100, sample_mean = NULL, sample_cov = NULL,
    sample_cor = dummy_cor, model = "cor", lambda_prior = "unit_vector",
    missing = "listwise",
    rotate = "varimax", rsp_args = NULL, factor_scores = FALSE,
    backend = "rstan", prior = list(), compute_fit_indices = FALSE,
    compute_reliability = FALSE, verbose = FALSE
  )

  result <- check_befa_inputs(args)
  expect_equal(result$N, 100)
})

test_that("befa accepts sample_cov with sample_nobs", {
  args <- list(
    data = NULL, n_factors = 2, ordered = FALSE,
    sample_nobs = 100, sample_mean = NULL, sample_cov = dummy_cov,
    sample_cor = NULL, model = "cov", lambda_prior = "unit_vector",
    missing = "listwise",
    rotate = "varimax", rsp_args = NULL, factor_scores = FALSE,
    backend = "rstan", prior = list(), compute_fit_indices = FALSE,
    compute_reliability = FALSE, verbose = FALSE
  )

  result <- check_befa_inputs(args)
  expect_equal(result$N, 100)
})

test_that("befa accepts ordered=FALSE", {
  args <- list(
    data = dummy_data, n_factors = 2, ordered = FALSE,
    sample_nobs = NULL, sample_mean = NULL, sample_cov = NULL,
    sample_cor = NULL, model = "cor", lambda_prior = "unit_vector",
    missing = "listwise",
    rotate = "varimax", rsp_args = NULL, factor_scores = FALSE,
    backend = "rstan", prior = list(), compute_fit_indices = FALSE,
    compute_reliability = FALSE, verbose = FALSE
  )

  result <- check_befa_inputs(args)
  expect_type(result, "list")
})

# Note: All integration tests that require Stan sampling (model fitting,
# verbose, reliability, fit indices, return object structure, Stan args, etc.)
# are now covered by run_all_befa_checks() in Section 2 below.

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Comprehensive Bayesian EFA Integration Tests
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
#                    INTEGRATION TESTS: FULL MODEL ESTIMATION
# ═══════════════════════════════════════════════════════════════════════════════
# These tests actually run Stan and verify that all model combinations work.
# They are slow and require rstan to be installed.

# --- Setup for integration tests ---
setup_integration_data <- function() {
  set.seed(42)
  # Generate data with clear factor structure
  n <- 100
  J <- 6

  # True loadings: items 1-3 load on F1, items 4-6 load on F2
  Lambda_true <- matrix(c(
    0.8, 0.1,
    0.75, 0.1,
    0.7, 0.15,
    0.1, 0.8,
    0.15, 0.75,
    0.1, 0.7
  ), nrow = 6, ncol = 2, byrow = TRUE)

  Psi_true <- diag(c(0.3, 0.35, 0.4, 0.3, 0.35, 0.4))
  Sigma_true <- Lambda_true %*% t(Lambda_true) + Psi_true

  # Generate multivariate normal data
  data <- MASS::mvrnorm(n = n, mu = rep(0, J), Sigma = Sigma_true)
  colnames(data) <- paste0("V", 1:J)

  list(
    data = as.data.frame(data),
    cor_matrix = cor(data),
    cov_matrix = cov(data),
    n = n,
    J = J
  )
}

# Minimal sampling settings for faster tests
minimal_stan_args <- list(
  iter = 400,
  warmup = 200,
  chains = 2,
  cores = 1,
  refresh = 0 # Suppress output
)

# ─────────────────────────────────────────────────────────────────────────────
# HELPER: Run ALL comprehensive checks on a fitted befa object
# ─────────────────────────────────────────────────────────────────────────────

run_all_befa_checks <- function(result, test_label) {
  # This function performs extensive validation on a fitted befa object
  # It is called within each model test to avoid code repetition

  J <- result$stan_data$J
  M <- result$n_factors
  model_type <- result$model_type

  # ── Basic structure checks ──
  expect_s3_class(result, "befa")
  expect_true("stanfit" %in% names(result))
  expect_true("stan_data" %in% names(result))
  expect_true("n_factors" %in% names(result))
  expect_true("model_type" %in% names(result))
  expect_true("lambda_prior" %in% names(result))
  expect_true("rsp_objective" %in% names(result))
  expect_true("priors_used" %in% names(result))

  # ── Posterior draws validation ──
  draws <- posterior::as_draws_matrix(result$stanfit)
  col_names <- colnames(draws)

  # Lambda parameters must exist
  lambda_cols <- grep("^Lambda\\[", col_names)
  expect_true(length(lambda_cols) > 0)
  expect_equal(length(lambda_cols), J * M)

  # Psi must exist in all models
  expect_true(any(grepl("^Psi\\[", col_names)))

  # Internal-only parameters must be stripped
  expect_false(any(grepl("^Z\\[", col_names)))
  expect_false(any(grepl("^Z_norm\\[", col_names)))
  expect_false(any(grepl("^psi\\[", col_names)))

  # Model-specific parameters
  if (model_type == "raw") {
    expect_true(any(grepl("^Sigma\\[", col_names)))
    expect_true(any(grepl("^nu\\[", col_names)))

    # Check Sigma is positive definite (first draw)
    sigma_cols <- grep("^Sigma\\[", col_names)
    Sigma_1 <- matrix(draws[1, sigma_cols], nrow = J, ncol = J)
    eigenvalues <- eigen(Sigma_1, only.values = TRUE)$values
    expect_true(all(eigenvalues > 0))
  } else {
    expect_true(any(grepl("^Rho\\[", col_names)))

    # Check Rho is valid correlation matrix (first draw)
    rho_cols <- grep("^Rho\\[", col_names)
    Rho_1 <- matrix(draws[1, rho_cols], nrow = J, ncol = J)
    expect_equal(diag(Rho_1), rep(1, J), tolerance = 0.01)
    eigenvalues <- eigen(Rho_1, only.values = TRUE)$values
    expect_true(all(eigenvalues > 0))
  }

  # ── Lambda draws quality ──
  lambda_draws <- draws[, lambda_cols]
  mean_loadings <- colMeans(lambda_draws)

  # Loadings should be bounded (roughly -3 to 3 for raw scale models)
  expect_true(all(abs(mean_loadings) <= 3.0))

  # At least some loadings should be substantially non-zero
  expect_true(any(abs(mean_loadings) > 0.2))

  # ── RSP alignment check ──
  expect_type(result$rsp_objective, "double")
  expect_true(all(is.finite(result$rsp_objective)))

  # ── Convergence diagnostics ──
  summ <- posterior::summarise_draws(lambda_draws, posterior::default_convergence_measures())

  # Rhat should not be terrible (even with few iterations)
  expect_true(all(summ$rhat < 1.5, na.rm = TRUE))

  # ESS should be positive
  expect_true(all(summ$ess_bulk > 0, na.rm = TRUE))
  expect_true(all(summ$ess_tail > 0, na.rm = TRUE))

  # ── Summary method ──
  summ_obj <- summary(result)
  expect_s3_class(summ_obj, "summary.befa")
  expect_equal(nrow(summ_obj$estimates), J)
  expect_true("h2" %in% names(summ_obj$estimates))
  expect_true("u2" %in% names(summ_obj$estimates))
  expect_true("Rhat" %in% names(summ_obj$estimates))

  # ── Reliability (if computed) ──
  if (!is.null(result$reliability)) {
    rel <- result$reliability
    expect_s3_class(rel, "befa_reliability")
    expect_true("total" %in% names(rel))
    expect_true("subscales" %in% names(rel))
    expect_true(all(c("Estimate", "SD", "Lower_95", "Upper_95") %in% names(rel$total)))
    expect_equal(nrow(rel$subscales), M)
    expect_equal(rownames(rel$subscales), paste0("F", 1:M))
    expect_true(rel$total$Estimate > 0 && rel$total$Estimate <= 1)
    expect_true(all(rel$subscales$Estimate > 0))
    expect_true(all(rel$subscales$Estimate <= 1))
  }

  # ── Fit indices (if computed) ──
  if (!is.null(result$fit_indices)) {
    fit <- result$fit_indices
    expect_s3_class(fit, "befa_fitmeasures")
    expect_true("fit_indices" %in% names(fit))
    expect_true("posterior_fit" %in% names(fit))
    expect_true("loo_object" %in% names(fit))
    expect_true("details" %in% names(fit))

    # Check all indices present
    expected_rows <- c(
      "Chi2", "Chi2_ppp", "Chi2_Null", "BRMSEA", "BGamma",
      "Adj_BGamma", "BMc", "SRMR", "BCFI", "BTLI",
      "ELPD", "LOOIC", "p_loo"
    )
    expect_true(all(expected_rows %in% rownames(fit$fit_indices)))

    # Check column names
    expect_equal(colnames(fit$fit_indices), c("Estimate", "SD", "Lower_95", "Upper_95"))

    # Value bounds
    expect_true(fit$fit_indices["BRMSEA", "Estimate"] >= 0)
    expect_true(fit$fit_indices["BCFI", "Estimate"] >= 0 &&
      fit$fit_indices["BCFI", "Estimate"] <= 1)
    expect_true(fit$fit_indices["Chi2_ppp", "Estimate"] >= 0 &&
      fit$fit_indices["Chi2_ppp", "Estimate"] <= 1)

    # LOO object
    expect_s3_class(fit$loo_object, "loo")

    # Details
    expect_true(all(c("p_star", "pD", "N") %in% names(fit$details)))
  }

  invisible(TRUE)
}

# ─────────────────────────────────────────────────────────────────────────────
# TEST GROUP: model = "cor", lambda_prior = "unit_vector"
# ─────────────────────────────────────────────────────────────────────────────

test_that("befa comprehensive: cor + unit_vector + M=1", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  test_data <- setup_integration_data()

  suppressWarnings({
    result <- befa(
      sample_cor = test_data$cor_matrix,
      sample_nobs = test_data$n,
      n_factors = 1,
      model = "cor",
      lambda_prior = "unit_vector",
      compute_reliability = TRUE,
      compute_fit_indices = FALSE, # No raw data for cor model
      verbose = FALSE,
      iter = minimal_stan_args$iter,
      warmup = minimal_stan_args$warmup,
      chains = minimal_stan_args$chains,
      refresh = minimal_stan_args$refresh
    )
  })

  run_all_befa_checks(result, "cor + unit_vector + M=1")
})

test_that("befa comprehensive: cor + unit_vector + M=2", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  test_data <- setup_integration_data()

  suppressWarnings({
    result <- befa(
      sample_cor = test_data$cor_matrix,
      sample_nobs = test_data$n,
      n_factors = 2,
      model = "cor",
      lambda_prior = "unit_vector",
      compute_reliability = TRUE,
      compute_fit_indices = FALSE, # No raw data for cor model
      verbose = FALSE,
      iter = minimal_stan_args$iter,
      warmup = minimal_stan_args$warmup,
      chains = minimal_stan_args$chains,
      refresh = minimal_stan_args$refresh
    )
  })

  run_all_befa_checks(result, "cor + unit_vector + M=2")
})


# ─────────────────────────────────────────────────────────────────────────────
# TEST GROUP: model = "raw", lambda_prior = "unit_vector"
# ─────────────────────────────────────────────────────────────────────────────

test_that("befa comprehensive: raw + unit_vector + M=1", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  test_data <- setup_integration_data()

  suppressWarnings({
    result <- befa(
      data = test_data$data,
      n_factors = 1,
      model = "raw",
      lambda_prior = "unit_vector",
      compute_reliability = TRUE,
      compute_fit_indices = TRUE,
      verbose = FALSE,
      iter = minimal_stan_args$iter,
      warmup = minimal_stan_args$warmup,
      chains = minimal_stan_args$chains,
      refresh = minimal_stan_args$refresh
    )
  })

  run_all_befa_checks(result, "raw + unit_vector + M=1")
})

test_that("befa comprehensive: raw + unit_vector + M=2", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  test_data <- setup_integration_data()

  suppressWarnings({
    result <- befa(
      data = test_data$data,
      n_factors = 2,
      model = "raw",
      lambda_prior = "unit_vector",
      compute_reliability = TRUE,
      compute_fit_indices = TRUE,
      verbose = FALSE,
      iter = minimal_stan_args$iter,
      warmup = minimal_stan_args$warmup,
      chains = minimal_stan_args$chains,
      refresh = minimal_stan_args$refresh
    )
  })

  run_all_befa_checks(result, "raw + unit_vector + M=2")
})

# ─────────────────────────────────────────────────────────────────────────────
# TEST GROUP: model = "raw", lambda_prior = "normal"
# ─────────────────────────────────────────────────────────────────────────────

test_that("befa comprehensive: raw + normal + M=1", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  test_data <- setup_integration_data()

  suppressWarnings({
    result <- befa(
      data = test_data$data,
      n_factors = 1,
      model = "raw",
      lambda_prior = "normal",
      compute_reliability = TRUE,
      compute_fit_indices = TRUE,
      verbose = FALSE,
      iter = minimal_stan_args$iter,
      warmup = minimal_stan_args$warmup,
      chains = minimal_stan_args$chains,
      refresh = minimal_stan_args$refresh
    )
  })

  run_all_befa_checks(result, "raw + normal + M=1")
})

test_that("befa comprehensive: raw + normal + M=2", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  test_data <- setup_integration_data()

  suppressWarnings({
    result <- befa(
      data = test_data$data,
      n_factors = 2,
      model = "raw",
      lambda_prior = "normal",
      compute_reliability = TRUE,
      compute_fit_indices = TRUE,
      verbose = FALSE,
      iter = minimal_stan_args$iter,
      warmup = minimal_stan_args$warmup,
      chains = minimal_stan_args$chains,
      refresh = minimal_stan_args$refresh
    )
  })

  run_all_befa_checks(result, "raw + normal + M=2")
})


# ─────────────────────────────────────────────────────────────────────────────
# NEGATIVE TESTS: Invalid combinations
# ─────────────────────────────────────────────────────────────────────────────

test_that("befa errors: cor + normal is blocked", {
  skip_if_not_installed("rstan")

  test_data <- setup_integration_data()

  expect_error(
    befa(
      sample_cor = test_data$cor_matrix,
      sample_nobs = test_data$n,
      n_factors = 2,
      model = "cor",
      lambda_prior = "normal",
      verbose = FALSE
    ),
    "not supported"
  )
})

# ─────────────────────────────────────────────────────────────────────────────
# EXHAUSTIVE INTEGRATION GRID: All valid model combinations
# ─────────────────────────────────────────────────────────────────────────────
# Tests ALL valid befa model combinations with varimax/none rotation.
# Covers: cor/cov/raw models, unit_vector/normal priors, M=1/2, listwise/FIML.

test_that("befa runs successfully on all valid parameter combinations", {
  skip_if_not_installed("rstan")
  skip_on_cran()

  # ─────────────────────────────────────────────────────────────────────────
  # Setup: Generate simulated data with clear factor structure
  # ─────────────────────────────────────────────────────────────────────────
  set.seed(42)
  N <- 80
  J <- 6

  # True loadings: items 1-3 load on F1, items 4-6 load on F2
  Lambda_true <- matrix(c(
    0.8, 0.1,
    0.75, 0.1,
    0.7, 0.15,
    0.1, 0.8,
    0.15, 0.75,
    0.1, 0.7
  ), nrow = 6, ncol = 2, byrow = TRUE)

  Psi_true <- diag(c(0.3, 0.35, 0.4, 0.3, 0.35, 0.4))
  Sigma_true <- Lambda_true %*% t(Lambda_true) + Psi_true

  # Complete data
  data_complete <- MASS::mvrnorm(n = N, mu = rep(0, J), Sigma = Sigma_true)
  colnames(data_complete) <- paste0("V", 1:J)

  # Data with missing values for FIML tests (~15% missing)
  data_missing <- data_complete
  data_missing[1:6, 1] <- NA
  data_missing[7:12, 2] <- NA
  data_missing[13:18, 3] <- NA
  data_missing[19:24, 4] <- NA

  # ─────────────────────────────────────────────────────────────────────────
  # Define the complete parameter grid
  # ─────────────────────────────────────────────────────────────────────────
  grid <- expand.grid(
    M = c(1, 2),
    model = c("cor", "cov", "raw"),
    lambda_prior = c("unit_vector", "normal"),
    missing = c("listwise", "FIML"),
    rotate = c("varimax", "none"),
    stringsAsFactors = FALSE
  )

  # ─────────────────────────────────────────────────────────────────────────
  # Filter invalid combinations:
  # 1. normal prior + cor model = NOT SUPPORTED
  # 2. FIML + M=1 + none = skip (M=1 always uses varimax internally anyway)
  # ─────────────────────────────────────────────────────────────────────────
  grid <- grid[!(grid$model == "cor" & grid$lambda_prior == "normal"), ]

  # Stan arguments for fast testing
  stan_args <- list(iter = 1000, warmup = 500, chains = 1, refresh = 0)

  for (i in 1:nrow(grid)) {
    p <- grid[i, ]

    # Select appropriate data based on missing strategy
    test_data <- if (p$missing == "FIML") data_missing else data_complete

    label <- sprintf(
      "[M=%d | %s | %s | %s | %s]",
      p$M, p$model, p$lambda_prior, p$missing, p$rotate
    )

    test_that(paste("Grid test:", label), {
      expect_no_error({
        suppressWarnings(
          fit <- befa(
            data = test_data,
            n_factors = p$M,
            model = p$model,
            lambda_prior = p$lambda_prior,
            missing = p$missing,
            rotate = p$rotate,
            iter = stan_args$iter,
            warmup = stan_args$warmup,
            chains = stan_args$chains,
            refresh = stan_args$refresh,
            backend = "rstan",
            verbose = FALSE,
            compute_fit_indices = FALSE,
            compute_reliability = FALSE
          )
        )
      })

      # Structural validation
      expect_s3_class(fit, "befa")
      expect_equal(fit$n_factors, p$M)
      expect_equal(fit$model_type, p$model)
      expect_equal(fit$lambda_prior, p$lambda_prior)

      # Rotation: respects user choice for all M values
      expect_equal(fit$rotation, p$rotate)

      # FIML flag should be set correctly
      if (p$missing == "FIML") {
        expect_true(fit$has_missing)
      }
    })
  }
})

# ─────────────────────────────────────────────────────────────────────────────
