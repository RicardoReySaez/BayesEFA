#' Check if Model is Unidimensional Forced (Internal)
#'
#' Determines if the model configuration forces the unidimensional ("uni") Stan model.
#' This occurs when M=1 and using unit_vector prior.
#'
#' @param n_factors Integer. Number of latent factors.
#' @param lambda_prior Character. The lambda prior type.
#' @return Logical. TRUE if the uni model should be used.
#' @keywords internal
#' @noRd
is_uni_model <- function(n_factors, lambda_prior) {
  n_factors == 1 && lambda_prior == "unit_vector"
}

#' Prepare Data Structure for Missing Value Handling (Internal)
#'
#' Separates complete and incomplete observations. For complete cases,
#' computes sufficient statistics (mean, covariance, correlation).
#' For incomplete cases, replaces NAs with sentinel value and precomputes
#' observation indices for efficient FIML in Stan.
#'
#' @param data Matrix. Raw data (N x J).
#' @param missing Character. "listwise" or "FIML".
#' @param model Character. "cor", "cov", or "raw".
#' @param sentinel Numeric. Value to represent missing in Stan (default -999).
#' @return A list with prepared data for Stan.
#' @keywords internal
#' @noRd
prepare_missing_data <- function(data, missing = "listwise", model = "raw", sentinel = -999) {
  data <- as.matrix(data)
  N <- nrow(data)
  J <- ncol(data)

  # Listwise: current behavior (remove all rows with any NA)
  if (missing == "listwise") {
    data_clean <- na.omit(data)
    N_clean <- nrow(data_clean)

    # Warning if rows were removed
    if (N_clean < N) {
      warning(sprintf(
        "Missing values detected. Observations reduced from %d to %d.",
        N, N_clean
      ), call. = FALSE)
    }

    # Compute sufficient statistics based on model type
    if (model == "raw") {
      m_obs <- unname(as.vector(colMeans(data_clean)))
      S_obs <- unname(stats::cov(data_clean) * (N_clean - 1) / N_clean)
      R_obs <- unname(stats::cor(data_clean))
    } else if (model == "cov") {
      # Covariance model: no means, but has variances
      m_obs <- rep(0, J)
      S_obs <- unname(stats::cov(data_clean) * (N_clean - 1) / N_clean)
      R_obs <- unname(stats::cor(data_clean))
    } else { # model == "cor"
      # Correlation model: no means, no variances
      m_obs <- rep(0, J)
      S_obs <- diag(J)
      R_obs <- unname(stats::cor(data_clean))
    }

    # Ensure Y is a clean matrix for Stan
    Y_clean <- unname(data_clean)
    attr(Y_clean, "na.action") <- NULL

    return(list(
      N_complete = N_clean,
      m_obs = m_obs,
      S_obs = S_obs,
      R_obs = R_obs,
      Y = Y_clean,
      N_incomplete = 0L,
      Y_miss = matrix(0, nrow = 0, ncol = J),
      has_missing = FALSE,
      sentinel = sentinel,
      # FIML placeholders (not used for listwise)
      max_obs = as.integer(J),
      n_obs = integer(0),
      obs_idx = matrix(0L, nrow = 0, ncol = J)
    ))
  }

  # FIML: For cor model, standardize ALL data first (column-wise, ignoring NAs)
  if (model == "cor") {
    data <- apply(data, 2, scale)
  }

  # Separate complete and incomplete observations
  complete_rows <- complete.cases(data)

  Y_complete <- unname(data[complete_rows, , drop = FALSE])
  Y_incomplete <- unname(data[!complete_rows, , drop = FALSE])

  N_complete <- nrow(Y_complete)
  N_incomplete <- nrow(Y_incomplete)

  # Compute sufficient statistics from complete cases only
  if (N_complete > 0) {
    if (model == "raw") {
      m_obs <- unname(as.vector(colMeans(Y_complete)))
      S_obs <- unname(stats::cov(Y_complete) * (N_complete - 1) / N_complete)
      R_obs <- unname(stats::cor(Y_complete))
    } else if (model == "cov") {
      m_obs <- rep(0, J)
      S_obs <- unname(stats::cov(Y_complete) * (N_complete - 1) / N_complete)
      R_obs <- unname(stats::cor(Y_complete))
    } else { # model == "cor"
      # Standardized model: data already scaled, compute correlation
      m_obs <- rep(0, J)
      S_obs <- diag(J)
      R_obs <- unname(stats::cor(Y_complete))
    }
  } else {
    # Edge case: ALL observations have missing values
    # Stan won't use sufficient statistics, pass placeholder zeros
    m_obs <- rep(0, J)
    S_obs <- diag(J)
    R_obs <- diag(J)

    warning(
      "All observations have at least one missing value. ",
      "FIML will use individual-level likelihood for all cases (slower).",
      call. = FALSE
    )
  }

  # Replace NAs with sentinel for Stan and compute FIML indices
  if (N_incomplete > 0) {
    Y_miss <- Y_incomplete

    # Compute precomputed indices for FIML (multiple indexing in Stan)
    n_obs <- apply(Y_incomplete, 1, function(row) sum(!is.na(row)))
    max_obs <- max(n_obs)

    # Build index matrix: obs_idx[i, k] = k-th observed column index for person i
    obs_idx <- matrix(0L, nrow = N_incomplete, ncol = max_obs)
    for (i in seq_len(N_incomplete)) {
      obs_cols <- which(!is.na(Y_incomplete[i, ]))
      obs_idx[i, seq_along(obs_cols)] <- obs_cols
    }

    # Replace NAs with sentinel for Stan
    Y_miss[is.na(Y_miss)] <- sentinel
  } else {
    Y_miss <- matrix(0, nrow = 0, ncol = J)
    n_obs <- integer(0)
    max_obs <- as.integer(J)
    obs_idx <- matrix(0L, nrow = 0, ncol = J)
  }

  return(list(
    N_complete = N_complete,
    m_obs = m_obs,
    S_obs = S_obs,
    R_obs = R_obs,
    Y = Y_complete,
    N_incomplete = N_incomplete,
    Y_miss = Y_miss,
    has_missing = (N_incomplete > 0),
    sentinel = sentinel,
    # FIML precomputed indices
    max_obs = as.integer(max_obs),
    n_obs = as.integer(n_obs),
    obs_idx = obs_idx
  ))
}

#' Input Validation for befa (Internal)
#' @param args A list containing all input arguments from the parent environment.
#' @return A list formatted for Stan.
#' @noRd
check_befa_inputs <- function(args) {
  # Extract all arguments explicitly
  data <- args$data
  n_factors <- args$n_factors
  ordered <- args$ordered
  sample_nobs <- args$sample_nobs
  sample_mean <- args$sample_mean
  sample_cov <- args$sample_cov
  sample_cor <- args$sample_cor
  model <- args$model
  missing <- args$missing
  factor_scores <- args$factor_scores
  compute_fit_indices <- args$compute_fit_indices
  compute_reliability <- args$compute_reliability
  verbose <- args$verbose

  # Helper function for logical argument validation
  check_is_logical <- function(x, name) {
    if (!is.logical(x) || length(x) != 1) {
      stop(sprintf("Argument '%s' must be a single logical value (TRUE/FALSE).", name), call. = FALSE)
    }
  }

  # Check logical arguments
  check_is_logical(ordered, "ordered")
  if (ordered) {
    stop("Ordinal/Categorical EFA (ordered = TRUE) is not yet implemented. Please use ordered = FALSE.", call. = FALSE)
  }
  check_is_logical(factor_scores, "factor_scores")
  check_is_logical(compute_fit_indices, "compute_fit_indices")
  check_is_logical(compute_reliability, "compute_reliability")
  check_is_logical(verbose, "verbose")

  # Validate 'missing' argument (default to "listwise" for backwards compatibility)
  missing_method <- if (is.null(missing)) "listwise" else missing
  if (!(missing_method %in% c("listwise", "FIML"))) {
    stop("Argument 'missing' must be 'listwise' or 'FIML'.", call. = FALSE)
  }

  # Check if data or summary statistics are available
  has_data <- !is.null(data)
  has_summary <- !is.null(sample_cov) || !is.null(sample_cor) || !is.null(sample_nobs)

  if (has_data && has_summary) {
    stop("Ambiguous input: Provide EITHER 'data' OR summary statistics, but not both.", call. = FALSE)
  }
  if (!has_data && !has_summary) {
    stop("No input provided: You must specify either 'data' or summary statistics.", call. = FALSE)
  }

  # FIML requires raw data
  if (missing_method == "FIML" && !has_data) {
    stop("FIML requires raw data. Cannot be used with summary statistics.", call. = FALSE)
  }

  # Validate raw data
  if (has_data) {
    if (!is.data.frame(data) && !is.matrix(data)) {
      stop("Argument 'data' must be a data.frame or matrix.", call. = FALSE)
    }

    # Check numeric values (before removing NAs to catch non-numeric)
    data_mat <- as.matrix(data)
    if (!all(apply(data_mat, 2, function(x) all(is.numeric(x) | is.na(x))))) {
      stop("Argument 'data' must contain numeric values only.", call. = FALSE)
    }
  }

  # Validate summary statistics
  if (has_summary) {
    if (is.null(sample_nobs) || !is.numeric(sample_nobs) || length(sample_nobs) != 1) {
      stop("When using summary statistics, 'sample_nobs' must be a single integer.", call. = FALSE)
    }

    if (!is.null(sample_cov) && !is.null(sample_cor)) {
      stop("Provide EITHER 'sample_cov' OR 'sample_cor', not both.", call. = FALSE)
    }

    # Helper for matrix validation
    validate_matrix <- function(mat, name) {
      if (!is.matrix(mat) || !is.numeric(mat)) stop(sprintf("'%s' must be a numeric matrix.", name), call. = FALSE)
      if (nrow(mat) != ncol(mat)) stop(sprintf("'%s' must be square.", name), call. = FALSE)
      if (!isSymmetric(unname(mat))) stop(sprintf("'%s' must be symmetric.", name), call. = FALSE)
    }

    if (!is.null(sample_cov)) validate_matrix(sample_cov, "sample_cov")
    if (!is.null(sample_cor)) validate_matrix(sample_cor, "sample_cor")
  }

  # Model dimensionality check
  J_check <- if (has_data) ncol(data) else max(ncol(sample_cov), ncol(sample_cor))

  if (is.null(n_factors) || !is.numeric(n_factors) || n_factors < 1) {
    stop("'n_factors' must be a positive integer.", call. = FALSE)
  }
  if (n_factors >= J_check) {
    stop("Number of factors must be smaller than the number of items (variables).", call. = FALSE)
  }

  # Prepare data based on missing value strategy
  if (has_data) {
    # Use prepare_missing_data for structured handling
    miss_data <- prepare_missing_data(data, missing = missing_method, model = model)

    # Build Stan data list
    final_list <- list(
      N_complete = miss_data$N_complete,
      N_incomplete = miss_data$N_incomplete,
      N = miss_data$N_complete, # Backwards compatibility
      J = J_check,
      M = n_factors,
      Y = miss_data$Y,
      Y_original = as.matrix(data), # Full data WITH NAs for FIML fit measures
      Y_miss = miss_data$Y_miss,
      has_missing = miss_data$has_missing,
      sentinel = miss_data$sentinel,
      # Always include all sufficient statistics (unified Stan model needs them)
      m_obs = miss_data$m_obs,
      S_obs = miss_data$S_obs,
      R_obs = miss_data$R_obs
    )

    # Add FIML precomputed indices if needed
    if (miss_data$has_missing && missing_method == "FIML") {
      final_list$max_obs <- miss_data$max_obs
      final_list$n_obs <- miss_data$n_obs
      final_list$obs_idx <- miss_data$obs_idx
    } else {
      # Placeholders for Stan (required but unused)
      final_list$max_obs <- as.integer(J_check)
      final_list$n_obs <- integer(0)
      final_list$obs_idx <- matrix(0L, nrow = 0, ncol = J_check)
    }
  } else {
    # Summary statistics path (no FIML possible, already checked above)
    J_check <- max(ncol(sample_cov), ncol(sample_cor))

    # Always provide all sufficient statistics (dummies for unused)
    if (!is.null(sample_cov)) {
      m_obs_val <- if (is.null(sample_mean)) rep(0, J_check) else as.vector(sample_mean)
      S_obs_val <- unname(sample_cov * (sample_nobs - 1) / sample_nobs)
      R_obs_val <- unname(stats::cov2cor(sample_cov))
    } else {
      m_obs_val <- rep(0, J_check)
      S_obs_val <- diag(J_check) # dummy
      R_obs_val <- unname(sample_cor)
    }

    final_list <- list(
      N_complete = sample_nobs,
      N_incomplete = 0L,
      N = sample_nobs, # Backwards compatibility
      J = J_check,
      M = n_factors,
      Y = matrix(0, nrow = 0, ncol = J_check), # Empty matrix for Stan
      Y_miss = matrix(0, nrow = 0, ncol = J_check),
      has_missing = FALSE,
      sentinel = -999,
      m_obs = m_obs_val,
      S_obs = S_obs_val,
      R_obs = R_obs_val,
      # FIML placeholders (required by Stan models even if not used)
      max_obs = as.integer(J_check),
      n_obs = integer(0),
      obs_idx = matrix(0L, nrow = 0, ncol = J_check)
    )
  }

  return(final_list)
}

#' Validate RSP Arguments (Internal)
#'
#' Validates and fills defaults for the Rotation-Sign-Permutation algorithm.
#'
#' @param rsp_args List or NULL. User provided arguments.
#' @return A named list with 'max_iter' and 'threshold'.
#' @keywords internal
#' @noRd
validate_rsp_args <- function(rsp_args) {
  # Default settings
  defaults <- list(max_iter = 1000, threshold = 1e-6)

  # If NULL, use defaults
  if (is.null(rsp_args)) {
    return(defaults)
  }

  # If user provides it, check
  if (!is.list(rsp_args)) {
    stop("Argument 'rsp_args' must be a named list (e.g., list(max_iter=1000)).", call. = FALSE)
  }

  # If names are wrong, the're just ignored
  valid_names <- names(defaults)
  unknowns <- setdiff(names(rsp_args), valid_names)
  if (length(unknowns) > 0) {
    warning(sprintf("Unknown arguments in 'rsp_args' ignored: %s", paste(unknowns, collapse = ", ")), call. = FALSE)
  }

  # Merge defaults with user input
  out <- defaults
  if (!is.null(rsp_args$max_iter)) out$max_iter <- rsp_args$max_iter
  if (!is.null(rsp_args$threshold)) out$threshold <- rsp_args$threshold

  # Checking plausible values
  if (!is.numeric(out$max_iter) || length(out$max_iter) != 1 || out$max_iter < 1) {
    stop("'rsp_args$max_iter' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(out$threshold) || length(out$threshold) != 1 || out$threshold <= 0) {
    stop("'rsp_args$threshold' must be a positive numeric value.", call. = FALSE)
  }

  return(out)
}

#' Validate Backend Selection (Internal)
#'
#' Checks if the requested backend is valid, available, and handles fallbacks.
#'
#' @param backend Character string. The selected backend.
#' @return Character string. The final backend to be used ("cmdstanr" or "rstan").
#' @keywords internal
#' @noRd
validate_backend <- function(backend) {
  # Rstan basic backend
  if (backend == "rstan") {
    return("rstan")
  }

  # Cmdstanr backend:
  if (backend == "cmdstanr") {
    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
      # Only warn once per session to avoid spamming the user
      if (is.null(getOption("BayesEFA.cmdstanr_warned"))) {
        msg <- paste0(
          "Backend 'cmdstanr' was requested (or default) but is not installed.\n",
          "Backend changed to 'rstan'.\n",
          "To install CmdStanR, visit: https://mc-stan.org/cmdstanr/articles/cmdstanr.html"
        )
        warning(msg, call. = FALSE)
        options(BayesEFA.cmdstanr_warned = TRUE)
      }

      return("rstan")
    }

    return("cmdstanr")
  }
}

#' Select the appropriate Stan model based on user input
#' @keywords internal
#' @noRd
select_befa_model <- function(model, lambda_prior, n_factors, verbose) {
  # Check if unidimensional model is forced
  is_uni <- is_uni_model(n_factors, lambda_prior)

  if (is_uni) {
    if (verbose) message("Note: For M=1 with UV, the specific 'uni' model (Beta-Scaled) is used.")
  }

  # "normal" prior is not supported for "cor" models
  if (lambda_prior == "normal" && model == "cor") {
    stop("The 'normal' prior is not supported for 'cor' models. Please use 'unit_vector'.", call. = FALSE)
  }

  # Unified model — flags are passed via stan_data
  return("befa_efa")
}

#' Validate Rotation Criterion (Internal)
#'
#' Validates that the rotation criterion is either 'varimax' or 'none'.
#'
#' @param rotate Character. The rotation criterion name.
#' @return Character. The validated rotation criterion ('varimax' or 'none').
#' @keywords internal
#' @noRd
validate_rotation_criterion <- function(rotate) {
  valid_rotations <- c("varimax", "none")
  rotate_clean <- tolower(trimws(rotate))

  if (!rotate_clean %in% valid_rotations) {
    stop(sprintf(
      "Invalid rotation '%s'. Valid options are: %s",
      rotate, paste(valid_rotations, collapse = ", ")
    ), call. = FALSE)
  }

  return(rotate_clean)
}


#' Validate loo_args (Internal)
#'
#' Validates and sets defaults for loo configuration.
#'
#' @param loo_args List of arguments.
#' @return List with validated arguments.
#' @keywords internal
#' @noRd
validate_loo_args <- function(loo_args) {
  # Defaults
  config <- list(
    r_eff = FALSE,
    cores = 1
  )

  if (is.null(loo_args)) {
    return(config)
  }
  if (!is.list(loo_args)) stop("Argument 'loo_args' must be a named list.", call. = FALSE)

  # Validate and update
  input_names <- names(loo_args)
  valid_names <- names(config)
  unknown_args <- setdiff(input_names, valid_names)
  if (length(unknown_args) > 0) {
    warning("Unknown arguments in 'loo_args' ignored: ", paste(unknown_args, collapse = ", "), call. = FALSE)
  }

  if (!is.null(loo_args$r_eff)) {
    if (!is.logical(loo_args$r_eff) || length(loo_args$r_eff) != 1) {
      stop("loo_args$r_eff must be a single logical value.", call. = FALSE)
    }
    config$r_eff <- loo_args$r_eff
  }

  if (!is.null(loo_args$cores)) {
    if (!is.numeric(loo_args$cores) || length(loo_args$cores) != 1 || loo_args$cores < 1) {
      stop("loo_args$cores must be a positive integer.", call. = FALSE)
    }
    config$cores <- as.integer(loo_args$cores)
  }

  return(config)
}
