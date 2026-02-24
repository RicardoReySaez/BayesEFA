#' Bayesian Fit Indices for Exploratory Factor Analysis
#'
#' Computes Bayesian analogs of common Structural Equation Modeling (SEM) fit
#' indices. This implementation follows the **DevM framework** proposed by
#' Garnier-Villarreal and Jorgensen (2020), which translates traditional
#' maximum-likelihood fit indices into the Bayesian framework. By evaluating
#' these indices across the MCMC draws, the function yields full posterior
#' distributions, allowing for the quantification of sampling uncertainty
#' (e.g., 95% credible intervals) rather than relying on a single point estimate.
#'
#' @param object A `befa` object returned by [befa()].
#' @param ... Additional arguments passed to Stan when automatically fitting the
#'   null model (e.g., `iter`, `chains`, `warmup`, `loo_args`).
#'
#' @details
#' **The DevM Framework (Deviance at the Posterior Mean)**
#'
#' In frequentist SEM, fit indices are derived from the exact fit test statistic
#' (\eqn{\chi^2}) and the model degrees of freedom (\eqn{df}). In the Bayesian
#' DevM framework (Garnier-Villarreal & Jorgensen, 2020), these quantities are
#' replaced by their Bayesian analogs:
#'
#' * **Observed Deviance (\eqn{D_i^{obs}})**: The model deviance evaluated at each
#'   posterior draw \eqn{i}, serving as the Bayesian analog of \eqn{\chi^2}.
#' * **Effective Number of Parameters (\eqn{pD})**: A penalty for model complexity.
#'   Following the authors' recommendations, this function estimates \eqn{pD}
#'   using the Leave-One-Out Information Criterion (\eqn{p_{loo}}) via Pareto Smoothed Importance
#' Sampling Leave-One-Out (PSIS-LOO, Vehtari et al., 2017).
#' * **Saturated Moments (\eqn{p^*})**: The number of nonredundant sample moments
#'   (e.g., \eqn{J \cdot (J + 3) / 2} for raw data, or \eqn{J \cdot (J + 1) / 2} for covariance matrices).
#' * **Bayesian Degrees of Freedom (\eqn{df_B})**: Defined as \eqn{df_B = p^* - pD}.
#' * **Bayesian Noncentrality Parameter (\eqn{\lambda_i})**: Defined at each iteration
#'   as \eqn{\lambda_i = D_i^{obs} - p^*}.
#'
#' **Absolute Fit Indices**
#'
#' Absolute indices evaluate how well the hypothesized model reproduces the
#' observed data. The function computes the following indices at each draw \eqn{i}:
#'
#' * **BRMSEA** (Bayesian Root Mean Square Error of Approximation):
#'   \deqn{\text{BRMSEA}_i = \sqrt{\max [0, (D_i^{obs} - p^*) / (df_B \cdot N)]}}
#' * **BGamma** (\eqn{\hat{\Gamma}}) and **Adjusted BGamma** (\eqn{\hat{\Gamma}_{adj}}):
#'   Bayesian analogs of the GFI.
#'   \deqn{\hat{\Gamma}_i = J / [J + (2 / N) \cdot (D_i^{obs} - p^*)]}
#'   \deqn{\hat{\Gamma}_{adj,i} = 1 - (p^* / df_B) \cdot (1 - \hat{\Gamma}_i)}
#' * **BMc** (Bayesian McDonald's Centrality Index):
#'   \deqn{\text{BMc}_i = \exp [-(1 / (2 \cdot N)) \cdot (D_i^{obs} - p^*)]}
#' * **Chi2_ppp** (Posterior Predictive P-value): The proportion of MCMC iterations
#'   where the deviance of replicated data exceeds the observed deviance
#'   \eqn{P(D^{rep} \geq D^{obs})}. Values closer to 0.50 indicate excellent fit,
#'   while values \eqn{< .05} suggest poor fit.
#' * **SRMR**: The Standardized Root Mean Square Residual, computed explicitly
#'   at each MCMC draw.
#'
#' **Incremental Fit Indices**
#'
#' Incremental indices compare the hypothesized model (\eqn{M_H}) against a worst-case
#' baseline independence model (\eqn{M_0}), where all observed variables are assumed
#' to be uncorrelated. *Note: The function automatically estimates this null model
#' under the hood using the same Stan settings, which may increase computation time.*
#'
#' * **BCFI** (Bayesian Comparative Fit Index):
#'   \deqn{\text{BCFI}_i = 1 - [(D_{H,i}^{obs} - p^*) / (D_{0,i}^{obs} - p^*)]}
#' * **BTLI** (Bayesian Tucker-Lewis Index / Non-Normed Fit Index):
#'   \deqn{\text{BTLI}_i = [ (D_{0,i}^{obs} - p^*) / df_{B0} - (D_{H,i}^{obs} - p^*) / df_{BH} ] / [ (D_{0,i}^{obs} - p^*) / df_{B0} - 1 ]}
#'
#' *(Note: Both BCFI and BTLI are strictly bounded to the \eqn{[0, 1]} interval in the final output).*
#'
#' **Predictive Fit Indices**
#'
#' Out-of-sample predictive accuracy is assessed using PSIS-LOO cross-validation:
#'
#' * **ELPD**: Expected Log Pointwise Predictive Density.
#' * **LOOIC**: Leave-One-Out Information Criterion (\eqn{-2 \cdot \text{ELPD}}).
#' * **p_loo**: Estimated effective number of parameters (\eqn{pD}).
#'
#' @return A list of class `befa_fitmeasures` containing:
#'   * `fit_draws`: A `draws_array` (iterations \eqn{\times} chains \eqn{\times} P)
#'     containing the posterior draws of all fit indices (`Chi2`, `Chi2_Null`,
#'     `BRMSEA`, `BGamma`, `Adj_BGamma`, `BMc`, `SRMR`, `BCFI`, `BTLI`).
#'   * `summary`: A `draws_summary` data frame computed via
#'     [posterior::summarise_draws()].
#'   * `posterior_fit`: A data frame containing the full MCMC posterior draws
#'     for all \eqn{\chi^2}-based indices, enabling custom plotting and density estimation.
#'   * `loo_object`: The full `loo` object returned by [loo::loo()].
#'   * `details`: A list with internal model constants (\eqn{p^*}, \eqn{pD}, \eqn{N},
#'     and `chi2_ppp`).
#'
#' @references
#' Garnier-Villarreal, M., & Jorgensen, T. D. (2020). Adapting Fit Indices for
#' Bayesian Structural Equation Modeling: Comparison to Maximum Likelihood.
#' *Psychological Methods, 25*(1), 46-70. <https://doi.org/10.1037/met0000224>
#'
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation
#' using leave-one-out cross-validation and WAIC. *Statistics and Computing, 27*(5),
#' 1413-1432. <https://doi.org/10.1007/s11222-016-9696-4>
#'
#' @examples
#' \dontrun{
#' # --------------------------------------------------------------- #
#' #    1. Bayesian Fit Measures after fitting Bayesian EFA model    #
#' # --------------------------------------------------------------- #
#'
#' # Fit Bayesian EFA model to the famous Grant-White School Data (Holzinger y Swineford , 1939)
#' befa_fit <- befa(
#'   data = HS_data,
#'   n_factors = 3,
#'   factor_scores = FALSE,
#'   compute_fit_indices = FALSE,
#'   compute_reliability = FALSE,
#'   backend = "cmdstanr",
#'   seed = 8
#' )
#'
#' # Compute Bayesian Fit Measures from a fitted befa object
#' bayesian_fit_measures <- befa_fit_measures(befa_fit)
#'
#' # Posterios summaries (NAs when posterior draws have sds close to zero)
#' bayesian_fit_measures$summary
#'
#' # # A tibble: 9 x 10
#' #   variable       mean   median      sd     mad       q5      q95  rhat ess_bulk ess_tail
#' #   <chr>         <dbl>    <dbl>   <dbl>   <dbl>    <dbl>    <dbl> <dbl>    <dbl>    <dbl>
#' # 1 Chi2        47.1     46.5    6.98    6.81     36.8     59.5     1.00    1410.    2596.
#' # 2 Chi2_Null  919.     919.     0       0       919.     919.     NA         NA       NA
#' # 3 BRMSEA       0.0543   0.0563 0.0204  0.0184    0.0159   0.0842  1.00    1418.    2596.
#' # 4 BGamma       0.992    0.992  0.00498 0.00495   0.983    0.999   1.00    1418.    2596.
#' # 5 Adj_BGamma   0.973    0.975  0.0163  0.0162    0.944    0.998   1.00    1418.    2596.
#' # 6 BMc          0.982    0.983  0.0112  0.0111    0.962    0.999   1.00    1418.    2596.
#' # 7 SRMR         0.0462   0.0454 0.00809 0.00777   0.0348   0.0609  1.00    1775.    2856.
#' # 8 BCFI         0.987    0.988  0.00779 0.00771   0.973    0.999   1.00    1418.    2596.
#' # 9 BTLI         0.989    1      0.0174  0         0.952    1       1.00    1587.      NA
#'
#' # Posterior draws of fit measures (transform to matrix)
#' posterior::as_draws_matrix(bayesian_fit_measures$fit_draws)
#'
#' # # A draws_matrix: 1000 iterations, 4 chains, and 9 variables
#' #     variable
#' # draw Chi2 Chi2_Null BRMSEA BGamma Adj_BGamma  BMc  SRMR BCFI
#' #   1    47       919  0.057   0.99       0.97 0.98 0.041 0.99
#' #   2    59       919  0.084   0.98       0.94 0.96 0.047 0.97
#' #   3    52       919  0.069   0.99       0.96 0.97 0.043 0.98
#' #   4    52       919  0.070   0.99       0.96 0.97 0.049 0.98
#' #   5    52       919  0.069   0.99       0.96 0.97 0.042 0.98
#' #   6    52       919  0.069   0.99       0.96 0.97 0.046 0.98
#' #   7    47       919  0.057   0.99       0.97 0.98 0.041 0.99
#' #   8    37       919  0.014   1.00       1.00 1.00 0.039 1.00
#' #   9    42       919  0.041   1.00       0.99 0.99 0.041 0.99
#' #   10   41       919  0.040   1.00       0.99 0.99 0.044 0.99
#' # # ... with 3990 more draws, and 1 more variables
#'
#' # We can also plot the posterior distribution of fit measures
#' hist(posterior::as_draws_matrix(bayesian_fit_measures$fit_draws)[, "SRMR"])
#'
#' # PSIS-LOO Expected log-predictive density value, uncertainty, and pareto-k values
#' bayesian_fit_measures$loo_object
#'
#' # Computed from 4000 by 301 log-likelihood matrix.
#' #
#' #          Estimate   SE
#' # elpd_loo  -3416.3 42.5
#' # p_loo        25.0  1.8
#' # looic      6832.6 85.0
#' # ------
#' # MCSE of elpd_loo is 0.1.
#' # MCSE and ESS estimates assume independent draws (r_eff=1).
#' #
#' # All Pareto k estimates are good (k < 0.7).
#' # See help('pareto-k-diagnostic') for details.
#'
#' # -----------------------------------------------------------------------------
#'
#' # -------------------------------------------------------------- #
#' #    2. Bayesian Fit Measures when fitting Bayesian EFA model    #
#' # -------------------------------------------------------------- #
#'
#' # Fit Bayesian EFA model to the famous Grant-White School Data (Holzinger y Swineford , 1939)
#' befa_fit <- befa(
#'   data = HS_data,
#'   n_factors = 3,
#'   factor_scores = FALSE,
#'   compute_fit_indices = FALSE,
#'   compute_reliability = TRUE,
#'   backend = "cmdstanr",
#'   seed = 8
#' )
#'
#' # befa_fit_measures output is inside the "fit_indices" object
#' # befa_fit$fit_indices
#'
#' # -----------------------------------------------------------------------------
#' }
#'
#' @export
befa_fit_measures <- function(object, ...) {
  # Class validation
  if (!inherits(object, "befa")) {
    stop("Input 'object' must be of class 'befa'.", call. = FALSE)
  }

  # Return cached results if already computed
  if (!is.null(object$fit_indices)) {
    return(object$fit_indices)
  }

  # Check if raw data exists in the object
  if (is.null(object$stan_data$Y_original) && is.null(object$stan_data$Y)) {
    warning("Fit measures require raw data!
             Since only summary statistics were provided, these will not be computed.",
      call. = FALSE
    )
    return(NULL) # Returns NULL so result$fit_indices becomes NULL
  }

  # Use original data with NAs for FIML-proper log-likelihood
  data <- as.matrix(object$stan_data$Y_original)
  if (is.null(object$stan_data$Y_original)) {
    # Fallback for older befa objects without Y_original
    data <- as.matrix(object$stan_data$Y)
  }

  # Compute internal metrics for the estimated model
  metrics <- compute_posterior_metrics(object, data)

  # Psis-loo estimate with proper r_eff calculation
  # Extract loo_args from '...' or object options
  args_list <- list(...)
  if (!is.null(args_list$loo_args)) {
    loo_config <- args_list$loo_args
  } else if (!is.null(object$options$loo_args)) {
    loo_config <- object$options$loo_args
  } else {
    loo_config <- list() # Will be filled with defaults by validator
  }
  loo_config <- validate_loo_args(loo_config)

  # log_lik is now a 3D array (iter x chains x N)
  if (loo_config$r_eff) {
    r_eff <- loo::relative_eff(exp(metrics$log_lik), cores = loo_config$cores)
  } else {
    r_eff <- NA
  }

  loo_res <- loo::loo(metrics$log_lik, r_eff = r_eff, cores = loo_config$cores)
  est_p_loo <- loo_res$estimates["p_loo", "Estimate"]

  # Compute null model metrics
  user_dots <- list(...)
  fit_obj <- object$stanfit

  # Auto-detect iteration settings if not provided
  if (is.null(user_dots$iter)) user_dots$iter <- fit_obj@sim$iter
  if (is.null(user_dots$chains)) user_dots$chains <- fit_obj@sim$chains
  if (is.null(user_dots$warmup)) user_dots$warmup <- fit_obj@sim$warmup

  null_stan_args <- normalize_stan_args(
    backend = "rstan", model_type = object$model_type,
    lambda_prior = object$lambda_prior,
    stan_data = object$stan_data, user_dots = user_dots
  )

  # Estimate null model (only in raw/cov cases; cor is just deterministic)
  if (object$model_type %in% c("raw", "cov")) message("Estimating Null Model...")

  null_metrics <- suppressWarnings(suppressMessages(
    compute_null_metrics(
      data = data,
      model_type = object$model_type,
      ll_saturated = metrics$ll_saturated,
      stan_data = object$stan_data,
      stan_args_list = null_stan_args,
      lambda_prior = object$lambda_prior,
      loo_config = loo_config
    )
  ))

  # Compute fit measures
  N <- metrics$N
  J <- metrics$J

  # Saturated Model Degrees of Freedom (p_star)
  if (object$model_type == "raw") {
    p_star <- J * (J + 3) / 2 # means + covariances
  } else if (object$model_type == "cov") {
    p_star <- J * (J + 1) / 2 # covariances only (no means)
  } else { # "cor"
    p_star <- J * (J - 1) / 2 # correlations only
  }

  # Hypothesized Model Stats
  lambda_hat <- pmax(0, metrics$chisq - p_star)
  df_bayes_H <- max(1, p_star - est_p_loo)

  # Null Model Stats
  lambda_null <- null_metrics$lambda_vec
  df_bayes_0 <- null_metrics$df

  # Absolute Fit Indices
  brmsea <- sqrt(lambda_hat / (df_bayes_H * N))
  bgamma <- J / (J + (2 / N) * lambda_hat)
  adj_bgamma <- 1 - (p_star / df_bayes_H) * (1 - bgamma)
  bmc <- exp(-0.5 * lambda_hat / N)

  # Incremental Fit Indices: BCFI and BTLI
  bcfi <- 1 - (lambda_hat / lambda_null)
  bcfi <- pmin(1, pmax(0, bcfi))

  ratio_null <- lambda_null / df_bayes_0
  ratio_H <- lambda_hat / df_bayes_H

  btli <- (ratio_null - ratio_H) / (ratio_null - 1)
  btli <- pmin(1, pmax(0, btli)) # Optional constraint, TLI can technically go > 1

  # Full posterior distribution of all fit measures
  posterior_df <- data.frame(
    Chi2        = metrics$chisq,
    Chi2_Null   = null_metrics$chisq_vec,
    BRMSEA      = brmsea,
    BGamma      = bgamma,
    Adj_BGamma  = adj_bgamma,
    BMc         = bmc,
    SRMR        = metrics$srmr,
    BCFI        = bcfi,
    BTLI        = btli
  )

  # Get chain structure from the stanfit for proper 3D array
  post_array <- posterior::as_draws_array(object$stanfit)
  n_iter <- dim(post_array)[1]
  n_chains <- dim(post_array)[2]
  fit_param_names <- colnames(posterior_df)
  n_fit_params <- length(fit_param_names)

  # Build 3D array: iterations x chains x parameters
  fit_array <- array(NA_real_,
    dim = c(n_iter, n_chains, n_fit_params),
    dimnames = list(
      iteration = seq_len(n_iter),
      chain = seq_len(n_chains),
      variable = fit_param_names
    )
  )

  # Fill the array by splitting the flat vectors back into chains
  for (ch in seq_len(n_chains)) {
    row_start <- (ch - 1) * n_iter + 1
    row_end <- ch * n_iter
    for (p in seq_len(n_fit_params)) {
      fit_array[, ch, p] <- as.numeric(posterior_df[row_start:row_end, p])
    }
  }

  # Convert to posterior draws_array
  fit_draws <- posterior::as_draws_array(fit_array)

  # Compute chi2_ppp (scalar, not a distribution)
  chi2_ppp <- mean(metrics$chisq_rep >= metrics$chisq)

  out <- list(
    fit_draws = fit_draws,
    summary = posterior::summarise_draws(fit_draws),
    posterior_fit = posterior_df,
    loo_object = loo_res,
    details = list(
      p_star   = p_star,
      pD       = est_p_loo,
      N        = N,
      chi2_ppp = chi2_ppp
    )
  )

  class(out) <- "befa_fitmeasures"
  return(out)
}
