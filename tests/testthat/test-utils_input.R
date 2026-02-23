# --- Setup común ---
get_default_args <- function(...) {
  defaults <- list(
    data = NULL, n_factors = 1, ordered = FALSE, sample_nobs = NULL,
    sample_mean = NULL, sample_cov = NULL, sample_cor = NULL,
    model = "cor", lambda_prior = "unit_vector", missing = "listwise",
    rotate = "varimax", rsp_args = NULL, factor_scores = FALSE,
    backend = "rstan", prior = list(), compute_fit_indices = FALSE,
    compute_reliability = FALSE, verbose = FALSE
  )
  utils::modifyList(defaults, list(...))
}

set.seed(123)
dummy_data <- as.data.frame(matrix(rnorm(100), ncol = 5))
dummy_cov <- cov(dummy_data)
dummy_cor <- cor(dummy_data)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Validate befa inputs
# ─────────────────────────────────────────────────────────────────────────────
# Custom functions used in this script
#
#
# ─────────────────────────────────────────────────────────────────────────────

test_that("check_befa_inputs errors on ambiguous input (data + summary)", {
  args <- get_default_args(data = dummy_data, sample_cov = dummy_cov, sample_nobs = 100)
  expect_error(check_befa_inputs(args), "Ambiguous input")
})

test_that("check_befa_inputs errors when no input provided", {
  args <- get_default_args() # Neither data nor summary
  expect_error(check_befa_inputs(args), "No input provided")
})

test_that("check_befa_inputs validates logical flags", {
  args <- get_default_args(data = dummy_data, ordered = "not_logical")
  expect_error(check_befa_inputs(args), "must be a single logical value")

  args2 <- get_default_args(data = dummy_data, factor_scores = 123)
  expect_error(check_befa_inputs(args2), "must be a single logical value")
})

test_that("check_befa_inputs rejects non-matrix/data.frame data", {
  args <- get_default_args(data = list(a = 1, b = 2))
  expect_error(check_befa_inputs(args), "must be a data.frame or matrix")
})

test_that("check_befa_inputs rejects non-numeric data", {
  bad_data <- data.frame(x = c("a", "b"), y = c("c", "d"))
  args <- get_default_args(data = bad_data, n_factors = 1)
  expect_error(check_befa_inputs(args), "must contain numeric values only")
})

test_that("check_befa_inputs warns on NA values in data", {
  data_na <- dummy_data
  data_na[1, 1] <- NA
  args <- get_default_args(data = data_na, model = "raw")
  expect_warning(check_befa_inputs(args), "Missing values detected")
})

test_that("check_befa_inputs errors when sample_nobs missing for summary stats", {
  args <- get_default_args(sample_cov = dummy_cov)
  expect_error(check_befa_inputs(args), "sample_nobs")
})

test_that("check_befa_inputs errors on both sample_cov and sample_cor", {
  args <- get_default_args(sample_cov = dummy_cov, sample_cor = dummy_cor, sample_nobs = 100)
  expect_error(check_befa_inputs(args), "EITHER 'sample_cov' OR 'sample_cor'")
})

test_that("check_befa_inputs validates matrix properties (square, symmetric)", {
  non_square <- matrix(1:6, nrow = 2)
  args <- get_default_args(sample_cov = non_square, sample_nobs = 100)
  expect_error(check_befa_inputs(args), "must be square")

  non_sym <- matrix(c(1, 2, 3, 4), nrow = 2)
  args2 <- get_default_args(sample_cov = non_sym, sample_nobs = 100)
  expect_error(check_befa_inputs(args2), "must be symmetric")
})

test_that("check_befa_inputs errors when n_factors >= J", {
  args <- get_default_args(data = dummy_data, n_factors = 5) # J=5
  expect_error(check_befa_inputs(args), "smaller than the number of items")
})

test_that("check_befa_inputs returns correct structure for 'cor' model", {
  args <- get_default_args(sample_cor = dummy_cor, sample_nobs = 100, n_factors = 2, model = "cor")
  result <- check_befa_inputs(args)
  expect_type(result, "list")
  expect_true(all(c("N_complete", "N_incomplete", "J", "M", "R_obs") %in% names(result)))
  expect_equal(result$N_complete, 100)
  expect_equal(result$N_incomplete, 0)
  expect_equal(result$J, 5)
})

test_that("check_befa_inputs returns correct structure for 'raw' model", {
  args <- get_default_args(data = dummy_data, n_factors = 2, model = "raw")
  result <- check_befa_inputs(args)
  expect_type(result, "list")
  expect_true(all(c("N_complete", "N_incomplete", "J", "M", "m_obs", "S_obs", "Y") %in% names(result)))
})

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1b: Missing data validation
# ─────────────────────────────────────────────────────────────────────────────

test_that("check_befa_inputs validates missing argument", {
  args <- get_default_args(data = dummy_data, missing = "invalid")
  expect_error(check_befa_inputs(args), "must be 'listwise' or 'FIML'")
})

test_that("check_befa_inputs errors on FIML with summary statistics", {
  args <- get_default_args(sample_cov = dummy_cov, sample_nobs = 100, missing = "FIML")
  expect_error(check_befa_inputs(args), "FIML requires raw data")
})

test_that("prepare_missing_data separates complete/incomplete correctly", {
  data_na <- dummy_data
  data_na[1:3, 1] <- NA
  result <- prepare_missing_data(data_na, "FIML", "raw")

  expect_equal(result$N_complete, nrow(dummy_data) - 3)
  expect_equal(result$N_incomplete, 3)
  expect_true(all(result$Y_miss[, 1] == -999))
  expect_true(result$has_missing)
})

test_that("prepare_missing_data computes sufficient stats from complete cases", {
  data_na <- dummy_data
  data_na[1, 1] <- NA
  result <- prepare_missing_data(data_na, "FIML", "raw")

  # Stats should match complete data (use unname for comparison)
  data_complete <- dummy_data[-1, ]
  expect_equal(result$m_obs, as.vector(colMeans(data_complete)))
})

test_that("prepare_missing_data listwise removes rows correctly", {
  data_na <- dummy_data
  data_na[1:3, 1] <- NA
  expect_warning(
    result <- prepare_missing_data(data_na, "listwise", "raw"),
    "Missing values detected"
  )

  expect_equal(result$N_complete, nrow(dummy_data) - 3)
  expect_equal(result$N_incomplete, 0)
  expect_false(result$has_missing)
  expect_equal(nrow(result$Y), nrow(dummy_data) - 3)
})

test_that("prepare_missing_data FIML with cor model uses correlation matrix", {
  data_na <- dummy_data
  data_na[1, 1] <- NA
  result <- prepare_missing_data(data_na, "FIML", "cor")

  # R_obs should be correlation matrix (diagonal = 1)
  expect_equal(diag(result$R_obs), rep(1, ncol(dummy_data)), tolerance = 1e-10)
})

test_that("prepare_missing_data handles all-missing edge case", {
  # Create data where every row has at least one NA
  data_all_missing <- dummy_data
  for (i in seq_len(nrow(data_all_missing))) {
    data_all_missing[i, (i %% 5) + 1] <- NA
  }

  expect_warning(
    result <- prepare_missing_data(data_all_missing, "FIML", "raw"),
    "All observations have at least one missing value"
  )

  expect_equal(result$N_complete, 0)
  expect_equal(result$N_incomplete, nrow(dummy_data))
  # Should return identity matrices as placeholders
  expect_equal(result$S_obs, diag(ncol(dummy_data)))
  expect_equal(result$R_obs, diag(ncol(dummy_data)))
})

test_that("prepare_missing_data replaces NAs with sentinel correctly", {
  data_na <- dummy_data
  data_na[1, 1] <- NA
  data_na[2, 3] <- NA
  result <- prepare_missing_data(data_na, "FIML", "raw", sentinel = -999)

  expect_equal(as.numeric(result$Y_miss[1, 1]), -999)
  expect_equal(as.numeric(result$Y_miss[2, 3]), -999)
  # Other values should NOT be sentinel
  expect_true(all(result$Y_miss[1, 2:5] != -999))
})

test_that("prepare_missing_data returns correct structure for FIML", {
  data_na <- dummy_data
  data_na[1, 1] <- NA
  result <- prepare_missing_data(data_na, "FIML", "raw")

  expected_names <- c(
    "N_complete", "m_obs", "S_obs", "R_obs", "Y",
    "N_incomplete", "Y_miss", "has_missing", "sentinel"
  )
  expect_true(all(expected_names %in% names(result)))
  expect_equal(result$sentinel, -999)
})

test_that("prepare_missing_data works with no missing values", {
  result <- prepare_missing_data(dummy_data, "FIML", "raw")

  expect_equal(result$N_complete, nrow(dummy_data))
  expect_equal(result$N_incomplete, 0)
  expect_false(result$has_missing)
  expect_equal(nrow(result$Y_miss), 0)
})

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Validate rsp arguments
# ─────────────────────────────────────────────────────────────────────────────

test_that("validate_rsp_args returns defaults when NULL", {
  result <- validate_rsp_args(NULL)
  expect_equal(result$max_iter, 1000)
  expect_equal(result$threshold, 1e-6)
})

test_that("validate_rsp_args merges user values with defaults", {
  result <- validate_rsp_args(list(max_iter = 500))
  expect_equal(result$max_iter, 500)
  expect_equal(result$threshold, 1e-6) # Default preserved
})

test_that("validate_rsp_args errors on non-list input", {
  expect_error(validate_rsp_args("not_a_list"), "must be a named list")
})

test_that("validate_rsp_args warns on unknown arguments", {
  expect_warning(validate_rsp_args(list(unknown_arg = 999)), "Unknown arguments")
})

test_that("validate_rsp_args errors on invalid max_iter", {
  expect_error(validate_rsp_args(list(max_iter = -1)), "positive integer")
  expect_error(validate_rsp_args(list(max_iter = "abc")), "positive integer")
})

test_that("validate_rsp_args errors on invalid threshold", {
  expect_error(validate_rsp_args(list(threshold = 0)), "positive numeric")
  expect_error(validate_rsp_args(list(threshold = -0.001)), "positive numeric")
})

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Validate Backend
# ─────────────────────────────────────────────────────────────────────────────

test_that("validate_backend returns 'rstan' when requested", {
  expect_equal(validate_backend("rstan"), "rstan")
})

test_that("validate_backend returns 'cmdstanr' when available", {
  skip_if_not_installed("cmdstanr")
  expect_equal(validate_backend("cmdstanr"), "cmdstanr")
})

test_that("validate_backend falls back to 'rstan' with warning if cmdstanr unavailable", {
  # This test may be tricky to run if cmdstanr IS installed.
  # Mock or skip based on environment.
  skip_if(requireNamespace("cmdstanr", quietly = TRUE), "cmdstanr is installed")
  expect_warning(result <- validate_backend("cmdstanr"), "not installed")
  expect_equal(result, "rstan")
})

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Validate befa model
# ─────────────────────────────────────────────────────────────────────────────

test_that("select_befa_model returns 'uni' suffix for M=1 with unit_vector", {
  result <- select_befa_model("cor", "unit_vector", n_factors = 1, verbose = FALSE)
  expect_equal(result, "befa_efa")
})

test_that("select_befa_model returns 'uv' for unit_vector with M > 1", {
  result <- select_befa_model("cor", "unit_vector", n_factors = 3, verbose = FALSE)
  expect_equal(result, "befa_efa")
})

test_that("select_befa_model returns 'norm' for normal prior", {
  result <- select_befa_model("raw", "normal", n_factors = 2, verbose = FALSE)
  expect_equal(result, "befa_efa")
})

test_that("select_befa_model blocks 'normal' prior with 'cor' model", {
  expect_error(
    select_befa_model("cor", "normal", n_factors = 2, verbose = FALSE),
    "not supported"
  )
})

test_that("select_befa_model prints message when verbose=TRUE and M=1", {
  expect_message(
    select_befa_model("cor", "unit_vector", n_factors = 1, verbose = TRUE),
    "uni"
  )
})

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Y_original Data Consistency
# ─────────────────────────────────────────────────────────────────────────────
# Data naming convention:
#   - Y: Data for Stan (complete cases for sufficient statistics)
#   - Y_original: Original raw data WITH NAs (for FIML fit measures, oblique rotations)
#   - Y_miss: Incomplete cases with sentinel for Stan FIML
# ─────────────────────────────────────────────────────────────────────────────

test_that("check_befa_inputs stores Y_original with complete data", {
  set.seed(123)
  data_complete <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)

  args <- get_default_args(data = data_complete, n_factors = 2, model = "raw")
  result <- check_befa_inputs(args)

  # Y_original should be identical to input
  expect_equal(nrow(result$Y_original), 100)
  expect_equal(result$Y_original, data_complete)
  # Y should also have all 100 rows (no NAs removed)
  expect_equal(nrow(result$Y), 100)
})


test_that("check_befa_inputs preserves NAs in Y_original with listwise", {
  set.seed(456)
  data_with_na <- matrix(rnorm(50 * 4), nrow = 50, ncol = 4)
  na_positions <- sample(length(data_with_na), 5)
  data_with_na[na_positions] <- NA

  args <- get_default_args(data = data_with_na, n_factors = 2, model = "raw")
  result <- suppressWarnings(check_befa_inputs(args))

  # Y_original preserves NAs
  expect_equal(nrow(result$Y_original), 50)
  expect_true(anyNA(result$Y_original))

  # Y has complete cases only
  n_complete <- sum(complete.cases(data_with_na))
  expect_equal(nrow(result$Y), n_complete)
  expect_false(anyNA(result$Y))
})


test_that("check_befa_inputs Y_original is NULL for summary statistics", {
  args <- get_default_args(
    sample_cov = diag(4),
    sample_nobs = 100,
    n_factors = 2,
    model = "raw"
  )

  result <- check_befa_inputs(args)

  # No raw data available
  expect_equal(nrow(result$Y), 0)
  expect_true(is.null(result$Y_original) || nrow(result$Y_original) == 0)
})

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Validate loo_args
# ─────────────────────────────────────────────────────────────────────────────

test_that("validate_loo_args returns defaults when NULL", {
  result <- validate_loo_args(NULL)
  expect_equal(result$r_eff, FALSE)
  expect_equal(result$cores, 1)
})

test_that("validate_loo_args updates defaults with user input", {
  result <- validate_loo_args(list(r_eff = FALSE, cores = 4))
  expect_equal(result$r_eff, FALSE)
  expect_equal(result$cores, 4)
})

test_that("validate_loo_args warns on unknown arguments", {
  expect_warning(validate_loo_args(list(unknown = "foo")), "Unknown arguments")
})

test_that("validate_loo_args errors on invalid input", {
  expect_error(validate_loo_args(list(r_eff = "yes")), "logical value")
  expect_error(validate_loo_args(list(cores = -1)), "positive integer")
})


# ─────────────────────────────────────────────────────────────────────────────
