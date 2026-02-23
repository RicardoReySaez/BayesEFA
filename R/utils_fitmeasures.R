#' Compute Log-Likelihood and Chi-Square (Internal)
#'
#' Extracts posterior 'Sigma' (or 'Rho') and 'nu' to compute pointwise log-likelihood
#' and Chi-square values relative to the saturated model.
#'
#' @param object A befa object.
#' @param data Matrix. Raw data (may contain NAs for FIML).
#' @return List with log_lik 3D array (iter x chains x N), chisq vector, and srmr vector.
#' @keywords internal
#' @noRd
compute_posterior_metrics <- function(object, data) {
  # Prepare objects
  N <- nrow(data)
  J <- ncol(data)
  model_type <- object$model_type
  has_missing <- anyNA(data)

  # Scale data for cor model
  if (model_type == "cor") {
    data <- apply(data, 2, scale)
  }

  # ---- FIML PRE-PROCESSING ----
  # Pre-split data into complete and incomplete subsets ONCE.
  # This drastically speeds up the log-likelihood loop by vectorizing complete cases.
  comp_idx <- which(stats::complete.cases(data))
  miss_idx <- which(!stats::complete.cases(data))

  data_fiml <- list(
    comp_idx = comp_idx,
    miss_idx = miss_idx,
    Y_comp   = data[comp_idx, , drop = FALSE],
    Y_miss   = data[miss_idx, , drop = FALSE]
  )

  # Saturated model: use EM if missing data, else direct computation
  if (has_missing) {
    # EM-based estimation for valid chi-square comparison with FIML
    h1_est <- estimate_h1_moments_em(data)
    M_sat <- h1_est$Mu
    S_sat <- h1_est$Sigma

    if (model_type == "cor") {
      # Standardize for correlation structure
      D_inv <- diag(1 / sqrt(diag(S_sat)))
      S_sat <- D_inv %*% S_sat %*% D_inv
      M_sat <- rep(0, J)
    } else if (model_type == "cov") {
      # Covariance model: set means to zero
      M_sat <- rep(0, J)
    }

    # Saturated LL using FIML (same N as hypothesized model)
    ll_saturated <- sum(compute_fiml_loglik(data_fiml, mu = M_sat, Sigma = S_sat))
  } else {
    # No missing: standard computation
    if (model_type == "cor") {
      S_sat <- stats::cor(data)
      M_sat <- rep(0, J)
    } else if (model_type == "cov") {
      S_sat <- stats::cov(data) * (N - 1) / N
      M_sat <- rep(0, J)
    } else { # raw
      S_sat <- stats::cov(data) * (N - 1) / N
      M_sat <- colMeans(data)
    }
    ll_saturated <- sum(mvnfast::dmvn(data, mu = M_sat, sigma = S_sat, log = TRUE))
  }

  # Extract posterior draws as 3D array (iter x chains x variables)
  draws_array <- posterior::as_draws_array(object$stanfit)
  n_iter <- dim(draws_array)[1]
  n_chains <- dim(draws_array)[2]
  n_draws <- n_iter * n_chains

  # Get variable names from the array
  var_names <- dimnames(draws_array)[[3]]

  # Identify columns based on model type
  if (model_type == "raw") {
    sig_idx <- grep("^Sigma\\[", var_names)
    nu_idx <- grep("^nu\\[", var_names)
    if (length(sig_idx) == 0) stop("Parameter 'Sigma' not found in Stan output.")
    if (length(nu_idx) == 0) stop("Parameter 'nu' not found in Stan output.")
  } else if (model_type == "cov") {
    sig_idx <- grep("^Sigma\\[", var_names)
    nu_idx <- NULL
    if (length(sig_idx) == 0) stop("Parameter 'Sigma' not found in Stan output.")
  } else { # cor
    sig_idx <- grep("^Rho\\[", var_names)
    nu_idx <- NULL
    if (length(sig_idx) == 0) stop("Parameter 'Rho' not found in Stan output.")
  }

  # Prepare 3D log-likelihood array (iter x chains x N) for loo::relative_eff()
  log_lik_array <- array(NA, dim = c(n_iter, n_chains, N))

  # Pre-allocate vectors for aggregate metrics (will be filled sequentially)
  chisq_vec <- chisq_rep_vec <- srmr_vec <- numeric(n_draws)

  # Pre-allocate loop variables
  Sigma_s <- matrix(0, J, J)
  Nu_s <- if (model_type == "raw") numeric(J) else rep(0, J) # cov and cor have means = 0

  # Main loop: iterate over chains and iterations separately
  draw_idx <- 0
  for (c in 1:n_chains) {
    for (i in 1:n_iter) {
      draw_idx <- draw_idx + 1

      # Reconstruct model-implied mean vector and cov/cor matrix
      Sigma_s <- matrix(draws_array[i, c, sig_idx], nrow = J, ncol = J)
      if (!is.null(nu_idx)) {
        Nu_s <- draws_array[i, c, nu_idx]
      }

      # Compute Pointwise Log-Likelihood (FIML-proper: handles NAs)
      ll_s <- compute_fiml_loglik(data_fiml, mu = Nu_s, Sigma = Sigma_s)
      log_lik_array[i, c, ] <- ll_s

      # Compute chi-square statistic: 2 * (LL_Saturated - LL_Model)
      chisq_vec[draw_idx] <- 2 * (ll_saturated - sum(ll_s))

      # Posterior predictive p-value (replicate N observations from model)
      Y_rep <- mvnfast::rmvn(n = N, mu = Nu_s, sigma = Sigma_s)
      ll_rep_model <- sum(mvnfast::dmvn(Y_rep, mu = Nu_s, sigma = Sigma_s, log = TRUE))

      # Observed moments for simulated data
      if (model_type == "cor") {
        S_sat_rep <- stats::cor(Y_rep)
        M_sat_rep <- rep(0, J)
      } else if (model_type == "cov") {
        S_sat_rep <- stats::cov(Y_rep) * (N - 1) / N
        M_sat_rep <- rep(0, J)
      } else { # raw
        S_sat_rep <- stats::cov(Y_rep) * (N - 1) / N
        M_sat_rep <- colMeans(Y_rep)
      }

      # Compute chi-square statistic
      ll_sat_rep <- sum(mvnfast::dmvn(Y_rep, mu = M_sat_rep, sigma = S_sat_rep, log = TRUE))
      chisq_rep_vec[draw_idx] <- 2 * (ll_sat_rep - ll_rep_model)

      # Prepare objects for SRMR statistic
      if (model_type == "raw" || model_type == "cov") {
        Imp_Cor <- stats::cov2cor(Sigma_s)
        Sample_Cor <- stats::cov2cor(S_sat)
      } else { # cor
        Imp_Cor <- Sigma_s
        Sample_Cor <- S_sat
      }

      # Compute SRMR
      diff_mat <- Sample_Cor - Imp_Cor
      resid_unique <- diff_mat[lower.tri(diff_mat, diag = TRUE)]
      srmr_vec[draw_idx] <- sqrt(mean(resid_unique^2))
    }
  }

  return(list(
    log_lik      = log_lik_array,
    chisq        = chisq_vec,
    chisq_rep    = chisq_rep_vec,
    srmr         = srmr_vec,
    ll_saturated = ll_saturated,
    N            = N,
    J            = J
  ))
}


#' Compute Null Model Metrics (Internal)
#'
#' Estimates the Null Model (Independence) to obtain the posterior distribution
#' of the Chi-square statistic and the baseline for incremental fit indices.
#'
#' @param data Matrix. Raw data (N x J).
#' @param model_type Character. "raw", "cov", or "cor".
#' @param ll_saturated Numeric. Log-likelihood of the saturated model.
#' @param stan_data List. The data list used for the main model.
#' @param ... Arguments passed to rstan::sampling (iter, chains, etc.).
#'
#' @param loo_config List. Configuration for LOO (r_eff, cores).
#'
#' @return A list with chisq vector, pD, df, and lambda vector.
#' @keywords internal
#' @noRd
compute_null_metrics <- function(data, model_type, ll_saturated, stan_data, stan_args_list, lambda_prior, loo_config = list(r_eff = TRUE, cores = 1)) {
  N <- nrow(data)
  J <- ncol(data)

  # Degrees of freedom based on model type
  if (model_type == "raw") {
    # means + covariances
    p_star <- J * (J + 3) / 2
  } else if (model_type == "cov") {
    # covariances only
    p_star <- J * (J + 1) / 2
  } else { # cor
    # correlations only
    p_star <- J * (J - 1) / 2
  }

  # Raw and cov null models require estimation via Stan
  if (model_type %in% c("raw", "cov")) {
    # Unified null model with model_type flag
    null_model <- stanmodels$befa_null

    # Prepare stan_data for null model: add model_type flag and required priors
    null_data <- stan_data
    null_data$model_type <- if (model_type == "raw") 1L else 2L

    # Always supply pr_nu and pr_sigma (unified model declares them both)
    if (is.null(null_data$pr_nu)) null_data$pr_nu <- c(0, 10)
    if (is.null(null_data$pr_sigma)) null_data$pr_sigma <- c(3, 0, 2.5)

    # Prepare sampling arguments
    call_args <- list(object = null_model)
    call_args$data <- null_data
    call_args <- c(call_args, stan_args_list)

    # Do not show anything by console!
    call_args$refresh <- 0
    call_args$show_messages <- FALSE
    call_args$open_progress <- FALSE

    # Sampling from the model
    fit_null <- do.call(rstan::sampling, call_args)

    # Posterior draws as 3D array (iter x chains x variables)
    draws_null <- posterior::as_draws_array(fit_null)
    n_iter <- dim(draws_null)[1]
    n_chains <- dim(draws_null)[2]
    n_draws <- n_iter * n_chains

    # Parameter indices
    var_names <- dimnames(draws_null)[[3]]
    sig_idx <- grep("^sigma\\[", var_names)

    # nu only exists in raw model
    if (model_type == "raw") {
      nu_idx <- grep("^nu\\[", var_names)
    }

    # Pre-allocate 3D log-likelihood array (iter x chains x N)
    log_lik_array <- array(NA, dim = c(n_iter, n_chains, N))
    chisq_vec <- numeric(n_draws)

    # ---- FIML PRE-PROCESSING ----
    comp_idx <- which(stats::complete.cases(data))
    miss_idx <- which(!stats::complete.cases(data))

    data_fiml <- list(
      comp_idx = comp_idx,
      miss_idx = miss_idx,
      Y_comp   = data[comp_idx, , drop = FALSE],
      Y_miss   = data[miss_idx, , drop = FALSE]
    )

    # Compute log-likelihood and chi-square statistic
    draw_idx <- 0
    for (c in 1:n_chains) {
      for (i in 1:n_iter) {
        draw_idx <- draw_idx + 1

        Sigma_s <- diag(draws_null[i, c, sig_idx]^2)

        if (model_type == "raw") {
          Nu_s <- draws_null[i, c, nu_idx]
        } else {
          Nu_s <- rep(0, J)
        }

        # Compute log-likelihood (FIML-proper: handles NAs)
        ll_s <- compute_fiml_loglik(data_fiml, mu = Nu_s, Sigma = Sigma_s)
        log_lik_array[i, c, ] <- ll_s

        # Compute chi-square statistic
        chisq_vec[draw_idx] <- 2 * (ll_saturated - sum(ll_s))
      }
    }

    # Compute r_eff and loo for null model
    if (!is.null(loo_config$r_eff) && loo_config$r_eff) {
      r_eff <- loo::relative_eff(exp(log_lik_array), cores = loo_config$cores)
    } else {
      r_eff <- NA
    }
    loo_null <- loo::loo(log_lik_array, r_eff = r_eff, cores = loo_config$cores)
    pD_null <- loo_null$estimates["p_loo", "Estimate"]

    # Final results
    df_bayes_0 <- max(1, p_star - pD_null)
    lambda_vec <- pmax(0, chisq_vec - p_star)

    return(list(
      chisq_vec  = chisq_vec,
      chisq_mean = mean(chisq_vec),
      pD         = pD_null,
      df         = df_bayes_0,
      lambda_vec = lambda_vec
    ))

    # Correlation model: no parameters (deterministic)
    # Null model assumes independence: each element ~ N(0,1)
    # Compute element-wise log-likelihood, removing NAs individually (not entire rows)
  } else { # cor
    # Standardize column-wise (handles NAs per column)
    z_data <- apply(data, 2, scale)

    # Element-wise log-likelihood: each z_ij ~ N(0, 1) independently
    ll_null <- sum(dnorm(na.omit(c(z_data)), mean = 0, sd = 1, log = TRUE))
    chisq_val <- 2 * (ll_saturated - ll_null)

    return(list(
      chisq_vec  = chisq_val,
      chisq_mean = chisq_val,
      pD         = 0,
      df         = max(1, p_star),
      lambda_vec = max(0, chisq_val - p_star)
    ))
  }
}

#' Print Bayesian Fit Measures
#'
#' Displays a formatted table of Bayesian fit indices from a `befa_fitmeasures` object.
#'
#' @param x A `befa_fitmeasures` object from [befa_fit_measures()].
#' @param digits Integer. Number of decimal places to display.
#' @param ... Ignored.
#'
#' @method print befa_fitmeasures
#' @export
print.befa_fitmeasures <- function(x, digits = 3, ...) {
  # --- Helper Formatting ---
  fmt <- function(n) {
    sapply(n, function(val) {
      if (is.na(val)) {
        return("")
      }
      sprintf(paste0("%.", digits, "f"), val)
    })
  }

  # Build and Prepare Data Frame
  fit_tab <- build_fit_summary_table(x)
  df_print <- as.data.frame(fit_tab)

  # Apply formatting to all columns
  df_print[] <- lapply(df_print, fmt)

  # Move Row names (Measure names) to the first column for the table printer
  df_print <- cbind(Index = rownames(df_print), df_print)
  rownames(df_print) <- NULL

  # Print Table
  cat("\n")
  # Capture the table width to adjust the notes afterwards
  table_width <- .print_styled_table(df_print, title = "Bayesian Fit Measures")

  # Construct Footer Notes
  notes <- c()

  # Standard Bayesian Notes
  notes <- c(notes, "Intervals represent 95% Highest Credible Intervals.")

  # LOO Specific Notes
  notes <- c(notes, "For ELPD/LOOIC, SD represents Standard Error.")
  notes <- c(notes, "p_loo and LOOIC derived from PSIS-LOO (loo package).")

  # Print Notes Wrapped to Table Width
  full_note <- paste("Note:", paste(notes, collapse = " "))
  cat(paste(strwrap(full_note, width = table_width), collapse = "\n"), "\n")

  invisible(x)
}
