#' Generate Default Priors (Internal)
#'
#' Creates the list of default hyperparameters based on the model combination.
#'
#' @param model Character. The data measurement model ("cor", "cov", or "raw").
#' @param lambda_prior Character. The factor loading prior ("unit_vector", "normal").
#' @param n_factors Integer. Number of factors.
#'
#' @return A named list with default values.
#' @keywords internal
#' @noRd
get_default_priors <- function(model, lambda_prior, n_factors) {
  defaults <- list()

  # Check if unidimensional model is forced
  is_uni <- is_uni_model(n_factors, lambda_prior)

  if (is_uni) {
    defaults$lambda <- c(1, 1) # Beta prior for uni model
  } else if (lambda_prior == "unit_vector") {
    defaults$xi <- 100
    defaults$h2 <- c(1, 1)
  } else if (lambda_prior == "normal") {
    defaults$lambda <- c(0, 1) # Normal prior for loadings
    defaults$psi <- c(0.5, 0.5)
  }

  # Model-specific priors
  if (model == "raw") {
    # Raw model: needs priors for means (nu) and residual SDs (sigma)
    defaults$nu <- c(0, 10)
    defaults$sigma <- c(3, 0, 2.5)
  } else if (model == "cov") {
    # Cov model: needs priors for residual SDs (sigma), but not means
    defaults$sigma <- c(3, 0, 2.5)
  }
  return(defaults)
}

#' Prepare and Validate Priors (Internal)
#'
#' Orchestrates the prior setup: Gets defaults, merges user input, validates final values.
#'
#' @param prior_user List. The 'prior' argument provided by the user.
#' @param model Character. "cor", "cov", or "raw".
#' @param lambda_prior Character.
#' @param n_factors Integer.
#'
#' @return A validated named list ready for Stan.
#' @keywords internal
#' @noRd
prepare_befa_priors <- function(prior_user, model, lambda_prior, n_factors) {
  # 1. Check if unidimensional model is forced
  is_uni <- is_uni_model(n_factors, lambda_prior)

  # 2. Generate default priors
  defaults <- get_default_priors(model, lambda_prior, n_factors)

  # 3. Check for valid priors and merge
  valid_names <- names(defaults)
  user_names <- names(prior_user)
  unknown_params <- setdiff(user_names, valid_names)

  if (length(unknown_params) > 0) {
    warning(sprintf(
      "The following parameters in 'prior' are not valid for this configuration and will be ignored: %s.",
      paste(unknown_params, collapse = ", ")
    ), call. = FALSE)
  }

  final_prior <- utils::modifyList(defaults, prior_user)

  # --- VALIDATIONS ---

  # Nu (only for RAW)
  if (model == "raw") {
    if (!is.numeric(final_prior$nu) || length(final_prior$nu) != 2) {
      stop("Prior 'nu' must be a vector of length 2.", call. = FALSE)
    }
    if (final_prior$nu[2] <= 0) stop("Prior 'nu': SD must be positive.", call. = FALSE)
  }

  # Sigma (for RAW and COV)
  if (model %in% c("raw", "cov")) {
    if (!is.numeric(final_prior$sigma) || length(final_prior$sigma) != 3) {
      stop("Prior 'sigma' must be a vector of length 3.", call. = FALSE)
    }
  }

  # EXCLUSIVE VALIDATION: UNI vs MULTIFACTORIAL
  if (is_uni) {
    # Validate lambda for the unidimensional model (Beta prior)
    if (!is.numeric(final_prior$lambda) || length(final_prior$lambda) != 2) {
      stop("Prior 'lambda' (uni model) must be a vector of length 2.", call. = FALSE)
    }
    if (any(final_prior$lambda <= 0)) {
      stop("Prior 'lambda' (uni model): Beta parameters must be positive.", call. = FALSE)
    }
  } else {
    # Multifactorial Validations
    if (lambda_prior == "unit_vector") {
      if (!is.numeric(final_prior$xi) || length(final_prior$xi) != 1) {
        stop("Prior 'xi' must be a single numeric value.", call. = FALSE)
      }
      if (!is.numeric(final_prior$h2) || length(final_prior$h2) != 2) {
        stop("Prior 'h2' must be a vector of length 2.", call. = FALSE)
      }
    }

    if (lambda_prior == "normal") {
      if (!is.numeric(final_prior$lambda) || length(final_prior$lambda) != 2) {
        stop("Prior 'lambda' must be a vector of length 2.", call. = FALSE)
      }
      if (!is.numeric(final_prior$psi) || length(final_prior$psi) != 2) {
        stop("Prior 'psi' must be a vector of length 2.", call. = FALSE)
      }
    }
  }

  return(final_prior)
}

#' Append Priors to Stan Data (Internal)
#'
#' Merges the data list with the prior hyperparameters, renaming them
#' to match the specific variable names expected by the Stan block.
#'
#' @param stan_data List. The output from check_befa_inputs (N, J, M, etc.).
#' @param final_prior List. The output from prepare_befa_priors.
#' @param model Character. "cor", "cov", or "raw".
#' @param lambda_prior Character. "unit_vector" or "normal".
#'
#' @return A single list containing both data and priors ready for sampling.
#' @keywords internal
#' @noRd
append_priors_to_data <- function(stan_data, final_prior, model, lambda_prior) {
  out <- stan_data
  n_factors <- stan_data$M
  is_uni <- is_uni_model(n_factors, lambda_prior)

  # ── Integer flags for the unified Stan model ──
  out$model_type <- switch(model,
    "raw" = 1L,
    "cov" = 2L,
    "cor" = 3L
  )

  if (is_uni) {
    out$lambda_type <- 2L
  } else {
    out$lambda_type <- switch(lambda_prior,
      "unit_vector" = 1L,
      "normal"      = 3L
    )
  }

  # nu prior (used by raw model only)
  out$pr_nu <- if (!is.null(final_prior$nu)) final_prior$nu else c(0, 10)

  # sigma prior (used by raw/cov models)
  out$pr_sigma <- if (!is.null(final_prior$sigma)) final_prior$sigma else c(3, 0, 2.5)

  # UV priors
  out$pr_xi <- if (!is.null(final_prior$xi)) final_prior$xi else 100
  out$pr_h2 <- if (!is.null(final_prior$h2)) final_prior$h2 else c(1, 1)

  # Lambda prior (shared by uni and normal; semantics differ)
  out$pr_Lambda <- if (!is.null(final_prior$lambda)) final_prior$lambda else c(0, 1)

  # Uniqueness prior (normal model only)
  out$pr_psi <- if (!is.null(final_prior$psi)) final_prior$psi else c(0.5, 0.5)

  return(out)
}
