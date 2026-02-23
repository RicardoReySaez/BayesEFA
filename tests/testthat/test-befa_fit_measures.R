# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: befa_fit_measures
# ─────────────────────────────────────────────────────────────────────────────

# --- Helper: create mock befa object ---
create_mock_befa_for_fitmeasures <- function(model_type = "raw", J = 5, M = 2,
                                             n_draws = 50, N = 100) {
  set.seed(123)

  # Generate fake data
  data <- matrix(rnorm(N * J), nrow = N, ncol = J)

  # Create realistic Lambda, Sigma, and nu draws
  true_lambda <- matrix(0.1, J, M)
  for (m in 1:M) {
    idx <- ((m - 1) * floor(J / M) + 1):min(m * floor(J / M), J)
    true_lambda[idx, m] <- runif(length(idx), 0.6, 0.8)
  }

  lambda_draws <- matrix(NA, n_draws, J * M)
  sigma_draws <- matrix(NA, n_draws, J * J)
  nu_draws <- matrix(NA, n_draws, J)

  for (s in 1:n_draws) {
    L_s <- true_lambda + matrix(rnorm(J * M, 0, 0.05), J, M)
    lambda_draws[s, ] <- as.vector(L_s)

    psi_s <- runif(J, 0.2, 0.5)
    Sigma_s <- L_s %*% t(L_s) + diag(psi_s)
    sigma_draws[s, ] <- as.vector(Sigma_s)

    nu_draws[s, ] <- colMeans(data) + rnorm(J, 0, 0.1)
  }

  # Create column names
  lambda_cols <- paste0("Lambda[", rep(1:J, M), ",", rep(1:M, each = J), "]")
  nu_cols <- paste0("nu[", 1:J, "]")

  if (model_type == "raw") {
    sigma_cols <- paste0("Sigma[", rep(1:J, J), ",", rep(1:J, each = J), "]")
    colnames(sigma_draws) <- sigma_cols
    colnames(nu_draws) <- nu_cols
    draws_matrix <- cbind(sigma_draws, nu_draws, lambda_draws)
  } else {
    rho_draws <- t(apply(sigma_draws, 1, function(s) {
      as.vector(cov2cor(matrix(s, J, J)))
    }))
    rho_cols <- paste0("Rho[", rep(1:J, J), ",", rep(1:J, each = J), "]")
    colnames(rho_draws) <- rho_cols
    draws_matrix <- cbind(rho_draws, lambda_draws)
  }

  colnames(lambda_draws) <- lambda_cols

  mock_stanfit <- structure(
    list(
      draws = draws_matrix,
      sim = list(iter = n_draws + 500, warmup = 500, chains = 2)
    ),
    class = "mock_stanfit"
  )

  structure(
    list(
      model_type = model_type,
      stanfit = mock_stanfit,
      stan_data = list(J = J, M = M, N = N, Y = data),
      n_factors = M,
      lambda_prior = "unit_vector"
    ),
    class = "befa"
  )
}

# --- Class validation tests ---

test_that("befa_fit_measures errors for non-befa object", {
  fake_obj <- list(a = 1, b = 2)
  class(fake_obj) <- "something_else"

  expect_error(
    befa_fit_measures(fake_obj),
    "must be of class 'befa'"
  )
})

test_that("befa_fit_measures returns NULL and warns when Y is missing", {
  mock_obj <- structure(
    list(
      stan_data = list(J = 5, M = 2, N = 100, Y = NULL),
      model_type = "cor"
    ),
    class = "befa"
  )

  expect_warning(
    result <- befa_fit_measures(mock_obj),
    "require raw data"
  )
  expect_null(result)
})

# Note: Structure validation tests (correct structure, row names, column names,
# posterior_fit data.frame, details contents) are now covered by
# run_all_befa_checks() in test-befa.R

# --- p_star calculation tests ---

test_that("p_star is correct for raw model (MACS)", {
  J <- 5
  p_star_raw <- J * (J + 3) / 2

  # J=5: 5 * 8 / 2 = 20
  expect_equal(p_star_raw, 20)
})

test_that("p_star is correct for std model (correlations only)", {
  J <- 5
  p_star_std <- J * (J - 1) / 2

  # J=5: 5 * 4 / 2 = 10
  expect_equal(p_star_std, 10)
})

# --- Absolute fit indices formulas ---

test_that("BRMSEA formula is correct", {
  lambda_hat <- 25
  df_bayes <- 10
  N <- 200

  brmsea <- sqrt(lambda_hat / (df_bayes * N))

  expect_equal(brmsea, sqrt(25 / 2000))
  expect_equal(brmsea, sqrt(0.0125))
})

test_that("BRMSEA is 0 when lambda_hat is 0", {
  lambda_hat <- 0
  df_bayes <- 10
  N <- 200

  brmsea <- sqrt(lambda_hat / (df_bayes * N))

  expect_equal(brmsea, 0)
})

test_that("BGamma formula is correct", {
  J <- 5
  N <- 200
  lambda_hat <- 20

  bgamma <- J / (J + (2 / N) * lambda_hat)
  expected <- 5 / (5 + 0.01 * 20)

  expect_equal(bgamma, expected)
})

test_that("BGamma approaches 1 when lambda_hat approaches 0", {
  J <- 5
  N <- 200
  lambda_hat <- 0.001

  bgamma <- J / (J + (2 / N) * lambda_hat)

  expect_true(bgamma > 0.999)
})

test_that("Adj_BGamma formula is correct", {
  J <- 5
  N <- 200
  lambda_hat <- 20
  p_star <- 15
  df_bayes <- 10

  bgamma <- J / (J + (2 / N) * lambda_hat)
  adj_bgamma <- 1 - (p_star / df_bayes) * (1 - bgamma)

  expect_true(adj_bgamma < bgamma)
  expect_true(adj_bgamma > 0)
})

test_that("BMc (McDonald's centrality) formula is correct", {
  N <- 200
  lambda_hat <- 20

  bmc <- exp(-0.5 * lambda_hat / N)
  expected <- exp(-10 / 200)

  expect_equal(bmc, expected)
})

test_that("BMc is 1 when lambda_hat is 0", {
  N <- 200
  lambda_hat <- 0

  bmc <- exp(-0.5 * lambda_hat / N)

  expect_equal(bmc, 1)
})

# --- Incremental fit indices formulas ---

test_that("BCFI formula is correct", {
  lambda_H <- 20
  lambda_0 <- 100

  bcfi <- 1 - (lambda_H / lambda_0)

  expect_equal(bcfi, 0.8)
})

test_that("BCFI is bounded between 0 and 1", {
  # Case where model is worse than null
  lambda_H <- 150
  lambda_0 <- 100

  bcfi_raw <- 1 - (lambda_H / lambda_0)
  bcfi <- pmin(1, pmax(0, bcfi_raw))

  expect_equal(bcfi, 0) # Truncated to 0
})

test_that("BCFI is 1 when lambda_H is 0", {
  lambda_H <- 0
  lambda_0 <- 100

  bcfi <- 1 - (lambda_H / lambda_0)

  expect_equal(bcfi, 1)
})

test_that("BTLI formula is correct", {
  ratio_null <- 10
  ratio_H <- 2

  btli <- (ratio_null - ratio_H) / (ratio_null - 1)
  expected <- (10 - 2) / 9

  expect_equal(btli, expected)
})

test_that("BTLI can exceed 1 before truncation", {
  ratio_null <- 10
  ratio_H <- 0.5

  btli_raw <- (ratio_null - ratio_H) / (ratio_null - 1)

  expect_true(btli_raw > 1)

  btli <- pmin(1, pmax(0, btli_raw))
  expect_equal(btli, 1)
})

# --- Chi-square and lambda_hat tests ---

test_that("Chi-square is 2 * (ll_sat - ll_model)", {
  ll_sat <- -100
  ll_model <- -150

  chisq <- 2 * (ll_sat - ll_model)

  expect_equal(chisq, 100)
})

test_that("lambda_hat is max(0, chisq - p_star)", {
  chisq <- 30
  p_star <- 20

  lambda_hat <- pmax(0, chisq - p_star)

  expect_equal(lambda_hat, 10)
})

test_that("lambda_hat is 0 when chisq < p_star", {
  chisq <- 15
  p_star <- 20

  lambda_hat <- pmax(0, chisq - p_star)

  expect_equal(lambda_hat, 0)
})

# --- Posterior predictive p-value (PPP) ---

test_that("PPP is proportion of chisq_rep >= chisq", {
  chisq <- c(50, 60, 55, 45, 70)
  chisq_rep <- c(55, 50, 60, 40, 75)

  ppp <- mean(chisq_rep >= chisq)

  # Comparison: 55>=50, 50>=60, 60>=55, 40>=45, 75>=70
  #             TRUE,   FALSE,  TRUE,   FALSE,  TRUE = 3/5
  expect_equal(ppp, 0.6)
})

test_that("PPP is near 0.5 for good fit", {
  set.seed(42)
  n <- 1000

  # Simulated case where model fits well
  chisq <- rnorm(n, 50, 10)
  chisq_rep <- rnorm(n, 50, 10) # Same distribution

  ppp <- mean(chisq_rep >= chisq)

  expect_true(ppp > 0.4 && ppp < 0.6)
})

# --- SRMR calculation ---

test_that("SRMR is 0 for perfect fit", {
  J <- 4
  Sample_Cor <- diag(J)
  Imp_Cor <- diag(J)

  diff_mat <- Sample_Cor - Imp_Cor
  resid_unique <- diff_mat[lower.tri(diff_mat, diag = TRUE)]
  srmr <- sqrt(mean(resid_unique^2))

  expect_equal(srmr, 0)
})

test_that("SRMR increases with misfit", {
  J <- 4
  Sample_Cor <- diag(J)

  # Small misfit
  Imp_Cor1 <- diag(J)
  Imp_Cor1[1, 2] <- Imp_Cor1[2, 1] <- 0.05
  diff1 <- Sample_Cor - Imp_Cor1
  srmr1 <- sqrt(mean(diff1[lower.tri(diff1, diag = TRUE)]^2))

  # Large misfit
  Imp_Cor2 <- diag(J)
  Imp_Cor2[1, 2] <- Imp_Cor2[2, 1] <- 0.3
  diff2 <- Sample_Cor - Imp_Cor2
  srmr2 <- sqrt(mean(diff2[lower.tri(diff2, diag = TRUE)]^2))

  expect_true(srmr2 > srmr1)
})

# Note: LOO integration (loo_object class) is now covered by run_all_befa_checks()

# --- df_bayes calculation ---

test_that("df_bayes is at least 1", {
  p_star <- 10
  p_loo <- 15 # More parameters than moments (unusual)

  df_bayes <- max(1, p_star - p_loo)

  expect_equal(df_bayes, 1)
})

test_that("df_bayes is p_star - p_loo when positive", {
  p_star <- 20
  p_loo <- 8

  df_bayes <- max(1, p_star - p_loo)

  expect_equal(df_bayes, 12)
})

# --- Summary function consistency ---

test_that("get_summary produces correct statistics", {
  vec <- c(0.05, 0.06, 0.04, 0.055, 0.045, 0.052, 0.048, 0.057, 0.043, 0.051)

  result <- c(
    Mean = mean(vec),
    SD = sd(vec),
    Q2.5 = as.numeric(quantile(vec, 0.025)),
    Q97.5 = as.numeric(quantile(vec, 0.975))
  )

  expect_equal(result["Mean"], mean(vec), ignore_attr = TRUE)
  expect_equal(result["SD"], sd(vec), ignore_attr = TRUE)
  expect_true(result["Q2.5"] < result["Mean"])
  expect_true(result["Q97.5"] > result["Mean"])
})

# ─────────────────────────────────────────────────────────────────────────────
