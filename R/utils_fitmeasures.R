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

  # Scale data for cor model or mean-center for cov model
  if (model_type == "cor") {
    data <- apply(data, 2, scale)
  } else if (model_type == "cov") {
    data <- scale(data, center = TRUE, scale = FALSE)
  }

  # ---- FIML PRE-PROCESSING ----
  # Pre-split data into complete and incomplete subsets ONCE.
  comp_idx <- which(stats::complete.cases(data))
  miss_idx <- which(!stats::complete.cases(data))

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
    sat_fiml <- list(
      comp_idx = comp_idx,
      miss_idx = miss_idx,
      Y_comp   = data[comp_idx, , drop = FALSE],
      Y_miss   = data[miss_idx, , drop = FALSE]
    )
    ll_saturated <- sum(compute_fiml_loglik(sat_fiml, mu = M_sat, Sigma = S_sat))
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

  # ── Prepare inputs for C++ backend ──

  # Convert model_type to integer flag (1 = raw, 2 = cov, 3 = cor)
  mt_int <- switch(model_type,
    "raw" = 1L,
    "cov" = 2L,
    "cor" = 3L
  )

  # 0-based indices for C++
  sig_idx_0 <- as.integer(sig_idx - 1L)
  nu_idx_0 <- if (!is.null(nu_idx)) as.integer(nu_idx - 1L) else integer(0)

  # FIML subsets for C++: 0-based indices
  comp_idx_0 <- as.integer(comp_idx - 1L)
  miss_idx_0 <- as.integer(miss_idx - 1L)

  # Observation patterns for missing rows (0-based)
  if (length(miss_idx) > 0) {
    obs_patterns_miss <- lapply(miss_idx, function(i) {
      as.integer(which(!is.na(data[i, ])) - 1L)
    })
  } else {
    obs_patterns_miss <- list()
  }

  # Complete / incomplete case subsets as matrices
  Y_comp <- if (length(comp_idx) > 0) data[comp_idx, , drop = FALSE] else matrix(0, nrow = 0, ncol = J)
  Y_miss <- if (length(miss_idx) > 0) data[miss_idx, , drop = FALSE] else matrix(0, nrow = 0, ncol = J)

  # Convert draws to a raw 3D array for C++
  draws_cube <- as.array(draws_array)

  # ── Call C++ backend ──
  cpp_result <- compute_posterior_metrics_cpp(
    post_draws        = draws_cube,
    data              = data,
    sig_idx           = sig_idx_0,
    nu_idx            = nu_idx_0,
    M_sat             = M_sat,
    S_sat             = S_sat,
    ll_saturated      = ll_saturated,
    model_type        = mt_int,
    has_missing       = has_missing,
    Y_comp            = Y_comp,
    Y_miss            = Y_miss,
    comp_idx_0        = comp_idx_0,
    miss_idx_0        = miss_idx_0,
    obs_patterns_miss = obs_patterns_miss
  )

  # Unpack results (match original return names)
  log_lik_array <- cpp_result$log_lik
  chisq_vec <- cpp_result$chisq
  chisq_rep_vec <- cpp_result$chisq_rep
  srmr_vec <- cpp_result$srmr

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
compute_null_metrics <- function(data, model_type, ll_saturated, stan_data, stan_args_list, lambda_prior, backend = "rstan", loo_config = list(r_eff = TRUE, cores = 1)) {
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
    # Prepare stan_data for null model: add model_type flag and required priors
    null_data <- stan_data
    null_data$model_type <- if (model_type == "raw") 1L else 2L

    # Always supply pr_nu and pr_sigma (unified model declares them both)
    if (is.null(null_data$pr_nu)) null_data$pr_nu <- c(0, 10)
    if (is.null(null_data$pr_sigma)) null_data$pr_sigma <- c(3, 0, 2.5)

    # Suppress messages when running the null model
    stan_args_list$refresh <- 0
    stan_args_list$show_message <- FALSE
    stan_args_list$show_messages <- FALSE
    
    # Exclude arguments we are about to explicitly pass to 'fit_befa_model'
    explicit_args <- c("model_name", "stan_data", "backend", "verbose", "model_type", "lambda_prior")
    stan_args_list <- stan_args_list[setdiff(names(stan_args_list), explicit_args)]

    # Sampling from the model
    # We use fit_befa_model to manage backend-specific behaviors securely
    fit_args <- c(
      list(
        model_name = "befa_null",
        stan_data = null_data,
        backend = backend,
        verbose = FALSE,
        model_type = model_type,
        lambda_prior = lambda_prior
      ),
      stan_args_list
    )
    fit_null <- suppressWarnings(suppressMessages(do.call(fit_befa_model, fit_args)))

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
    # Mean-center data for cov model before computing pointwise log-likelihood
    if (model_type == "cov") {
      data <- scale(data, center = TRUE, scale = FALSE)
    }

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

        Sigma_s <- diag(as.numeric(draws_null[i, c, sig_idx])^2)

        if (model_type == "raw") {
          Nu_s <- as.numeric(draws_null[i, c, nu_idx])
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
    loo_null <- suppressWarnings(
      loo::loo(log_lik_array, r_eff = r_eff, cores = loo_config$cores)
    )
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
