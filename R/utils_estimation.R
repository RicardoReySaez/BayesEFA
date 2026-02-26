#' Fit BEFA Model (Internal)
#'
#' Handles the backend selection, verbose messaging, and sampling execution.
#'
#' @param model_name Character. Name of the pre-compiled stanmodel (e.g., "befa_efa").
#' @param stan_data List. The data prepared for Stan.
#' @param backend Character. "cmdstanr" or "rstan".
#' @param verbose Logical. Whether to print start-up messages.
#' @param model_type Character. "cor", "cov", or "raw" (for printing).
#' @param lambda_prior Character. Prior name (for printing).
#' @param ... Additional arguments passed to sampling (iter, chains, etc.).
#'
#' @return A \code{stanfit} object.
#' @keywords internal
#' @noRd
fit_befa_model <- function(model_name, stan_data, backend, verbose, model_type, lambda_prior, ...) {
  # Check stan model existence
  if (!exists(model_name, where = stanmodels)) {
    stop(sprintf("Stan model '%s' not found internally.", model_name), call. = FALSE)
  }
  rstan_obj <- stanmodels[[model_name]]

  # Messages for users
  if (verbose) {
    txt_title <- "Bayesian Exploratory Factor Analysis (BEFA)"
    N_total <- stan_data$N_complete + stan_data$N_incomplete
    txt_dims <- sprintf("   Data  : N = %d observations | J = %d items", N_total, stan_data$J)
    txt_spec <- sprintf("   Model : %d Factors | %s | %s Prior", stan_data$M, toupper(model_type), toupper(lambda_prior))
    txt_eng <- sprintf("   Engine: %s", ifelse(backend == "cmdstanr", "CmdStanR", "RStan"))

    width <- max(nchar(txt_title), nchar(txt_dims), nchar(txt_spec), nchar(txt_eng)) + 2
    thick_line <- paste0(rep("\u2017", width), collapse = "")

    message("\n", txt_title)
    message(thick_line)
    message(txt_spec)
    message(txt_dims)
    message(txt_eng)
    message(thick_line)
    message("   Sampling posterior draws...\n")
  }

  # Normalize arguments
  stan_args <- normalize_stan_args(backend, model_type, lambda_prior, stan_data, list(...))
  stan_data_for_sampling <- stan_data
  stan_data_for_sampling$Y_original <- NULL
  stan_data_for_sampling$Y <- NULL

  stan_args$data <- stan_data_for_sampling

  out <- tryCatch(
    {
      # Cmdstanr estimation: save model on cache
      if (backend == "cmdstanr") {
        # Save .exe files in cache
        cache_dir <- tools::R_user_dir("BayesEFA", which = "cache")
        if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
        # Get Stan model code and write cmdstanr file
        stan_code <- rstan::get_stancode(rstan_obj)
        stan_file <- cmdstanr::write_stan_file(
          code = stan_code,
          dir = cache_dir,
          basename = paste0(model_name, "_hash")
        )
        # Compile cmdstanr .exe file
        mod <- cmdstanr::cmdstan_model(stan_file, quiet = TRUE)
        # Sampling from the model
        fit_cmd <- do.call(mod$sample, stan_args)
        # Make the output equivalent to rstan
        brms::read_csv_as_stanfit(fit_cmd$output_files())

        # Rstan branch: simple
      } else {
        stan_args$object <- rstan_obj
        do.call(rstan::sampling, stan_args)
      }
    },
    error = function(e) {
      stop("Stan sampling failed: ", e$message, call. = FALSE)
    }
  )

  return(out)
}

#' Normalize Stan Arguments (Internal)
#'
#' Harmonizes arguments between RStan and CmdStanR to ensure consistent behavior.
#'
#' @param backend Character. "rstan" or "cmdstanr".
#' @param user_dots List. The captured \code{...} from the main function.
#' @return A named list of arguments ready for sampling.
#' @keywords internal
#' @noRd
normalize_stan_args <- function(backend, model_type, lambda_prior, stan_data, user_dots) {
  # Package defaults for both backends
  defaults <- list(
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = 1,
    adapt_delta = 0.95,
    max_treedepth = 15
  )

  args <- user_dots

  # Make equivalent inputs between backends
  # Map CmdStanR terminology to RStan (Total Iterations)
  if (!is.null(args$iter_sampling) && !is.null(args$iter_warmup)) {
    args$warmup <- args$iter_warmup
    args$iter <- args$iter_sampling + args$iter_warmup
  }

  # Apply defaults or calculate missing values
  iter <- if (!is.null(args$iter)) args$iter else defaults$iter
  warmup <- if (!is.null(args$warmup)) args$warmup else defaults$warmup

  # Fallback: standard Stan behavior (warmup is 50% of iter)
  if (!is.null(args$iter) && is.null(args$warmup)) {
    warmup <- floor(iter / 2)
  }

  # Refresh values
  if (!is.null(args$refresh)) {
    refresh <- args$refresh
  } else {
    # Default: 25% of total iterations (updates 4 times per chain)
    refresh <- max(1, floor(iter / 4))
  }

  # Map Cores / Parallel Chains
  chains <- if (!is.null(args$chains)) args$chains else defaults$chains

  cores_input <- NULL
  if (!is.null(args$cores)) cores_input <- args$cores
  if (!is.null(args$parallel_chains)) cores_input <- args$parallel_chains

  # Non-parallelized estimation if not indicated
  cores <- if (!is.null(cores_input)) cores_input else 1

  # Control parameters (extract from list or direct args)
  adapt_delta <- defaults$adapt_delta
  max_treedepth <- defaults$max_treedepth

  if (!is.null(args$control)) {
    if (!is.null(args$control$adapt_delta)) adapt_delta <- args$control$adapt_delta
    if (!is.null(args$control$max_treedepth)) max_treedepth <- args$control$max_treedepth
  }
  if (!is.null(args$adapt_delta)) adapt_delta <- args$adapt_delta
  if (!is.null(args$max_treedepth)) max_treedepth <- args$max_treedepth

  # Ensure reproducibility of initial values if seed is provided
  if (!is.null(args$seed)) {
    withr::local_seed(args$seed)
  }

  if (is.null(args$init)) {
    inits <- get_befa_inits(model_type, lambda_prior, stan_data, chains)
  } else {
    inits <- user_dots$init
  }

  # Format Output by Backend
  final_args <- list()
  if (!is.null(args$seed)) final_args$seed <- args$seed

  # Identify known keys to isolate extra arguments
  known_keys <- c(
    "iter", "warmup", "chains", "cores", "control",
    "iter_sampling", "iter_warmup", "parallel_chains",
    "adapt_delta", "max_treedepth", "seed", "init", "refresh",
    "show_messages", "show_message", "show_exceptions",
    "loo_args"
  )
  extra_args <- args[!names(args) %in% known_keys]

  if (backend == "rstan") {
    # RStan requirements
    final_args$chains <- chains
    final_args$iter <- iter
    final_args$warmup <- warmup
    final_args$cores <- cores
    final_args$refresh <- refresh
    final_args$init <- inits

    # RStan expects a 'control' list
    ctrl_list <- if (!is.null(args$control)) args$control else list()
    ctrl_list$adapt_delta <- adapt_delta
    ctrl_list$max_treedepth <- max_treedepth
    final_args$control <- ctrl_list

    # Handle message visibility (RStan uses plural 'show_messages')
    if (!is.null(args$show_messages)) {
      final_args$show_messages <- args$show_messages
    } else if (!is.null(args$show_message)) {
      final_args$show_messages <- args$show_message
    }

    final_args <- c(final_args, extra_args)
  } else {
    # CmdStanR requirements
    final_args$chains <- chains
    final_args$iter_warmup <- warmup
    final_args$iter_sampling <- iter - warmup
    final_args$parallel_chains <- cores
    final_args$refresh <- refresh
    final_args$init <- inits

    # CmdStanR accepts direct arguments
    final_args$adapt_delta <- adapt_delta
    final_args$max_treedepth <- max_treedepth

    # Handle message/exception visibility
    # Map 'show_messages' (if present) to singular 'show_message' for cmdstanr
    if (!is.null(args$show_message)) {
      final_args$show_message <- args$show_message
    } else if (!is.null(args$show_messages)) {
      final_args$show_message <- args$show_messages
    }

    if (!is.null(args$show_exceptions)) {
      final_args$show_exceptions <- args$show_exceptions
    }

    final_args <- c(final_args, extra_args)
  }

  return(final_args)
}

#' Generate informative starting values for BEFA models
#'
#' @param model_type Character. "raw", "cov", "cor", or "null".
#' @param stan_data List. The data list sent to Stan.
#' @param n_chains Integer. Number of Markov chains.
#' @return A list of lists containing starting values.
#' @keywords internal
#' @noRd
get_befa_inits <- function(model_type, lambda_prior, stan_data, n_chains) {
  # If stan_data arrives malformed due to a flow error, this prevents an ugly crash
  if (!is.list(stan_data)) {
    return(NULL)
  }

  J <- stan_data$J
  M <- stan_data$M

  # Base statistics
  m_obs <- if (!is.null(stan_data[["m_obs"]])) stan_data[["m_obs"]] else rep(0, J)

  # Ensure S_obs is a matrix for diag()
  S_mat <- if (!is.null(stan_data[["S_obs"]])) as.matrix(stan_data[["S_obs"]]) else diag(J)
  s_obs <- sqrt(diag(S_mat))

  # Initial communalities (h2)
  R <- stats::cov2cor(S_mat)
  diag(R) <- 0
  initial_h2 <- pmin(0.9, pmax(0.1, apply(abs(R), 1, max)))

  # Determine lambda_type to match Stan logic
  # 1 = UV, 2 = uni, 3 = normal
  is_uni <- (M == 1 && lambda_prior == "unit_vector")
  lambda_type <- if (is_uni) 2 else if (lambda_prior == "unit_vector") 1 else 3

  lapply(1:n_chains, function(i) {
    inits <- list()

    # nu (Means): Only present if model_type == "raw" (1)
    if (model_type == "raw") {
      inits$nu <- as.array(m_obs + rnorm(J, 0, 0.1))
    }

    # sigma (Scale): Present if model_type <= 2 (raw or cov)
    if (model_type %in% c("raw", "cov")) {
      inits$sigma <- as.array(s_obs * exp(rnorm(J, 0, 0.1)))
    }

    # Factor Parameters (Z, h2, Lambda...)
    # For befa_null (null model), all factor params are zero-sized -> omit all.
    # For befa_efa, only provide inits for the active parameterization.
    if (model_type != "null") {
      # Type 1: Unit Vector (UV) -> h2, Z
      if (lambda_type == 1) {
        inits$h2 <- pmin(0.95, pmax(0.05, as.array(initial_h2 + rnorm(J, 0, 0.05))))
        inits$Z <- matrix(rnorm(J * M, 0, 1), nrow = J, ncol = M)

        # Type 2: Unidimensional -> Lambda_uni
      } else if (lambda_type == 2) {
        inits$Lambda_uni <- matrix(
          pmin(0.95, pmax(-0.95, sqrt(initial_h2) + rnorm(J, 0, 0.05))),
          nrow = J, ncol = 1
        )

        # Type 3: Normal -> Lambda_norm, psi
      } else {
        inits$Lambda_norm <- matrix(sqrt(initial_h2) + rnorm(J * M, 0, 0.1), nrow = J, ncol = M)
        inits$psi <- as.array((1 - initial_h2) * (s_obs^2) * exp(rnorm(J, 0, 0.1)))
      }
    }

    return(inits)
  })
}

#' Strip Internal Parameters from Stan Output (Internal)
#'
#' Removes parameters that are internal to the Stan model and not intended
#' for the end user. This keeps the stanfit object clean and consistent.
#'
#' For unit_vector models: removes `Z` and `Z_norm` (UV construction internals).
#' For normal models: removes `psi` (replaced by the `Psi` matrix from generated quantities).
#'
#' @param stanfit A stanfit object.
#' @param lambda_prior Character. "unit_vector" or "normal".
#' @param n_factors Integer. Number of factors (to determine if uni model).
#' @return The stanfit object with internal parameters removed.
#' @keywords internal
#' @noRd
strip_stan_params <- function(stanfit, lambda_prior, n_factors) {
  is_uni <- is_uni_model(n_factors, lambda_prior)

  # Define parameters to strip based on prior
  if (is_uni) {
    # Uni models have Lambda_uni — strip it (keep Lambda which is a copy)
    params_to_strip <- c("Lambda_uni")
  } else if (lambda_prior == "unit_vector") {
    params_to_strip <- c("Z", "Z_norm")
  } else if (lambda_prior == "normal") {
    # Strip Lambda_norm and psi (user gets Lambda and Psi)
    params_to_strip <- c("Lambda_norm", "psi")
  } else {
    return(stanfit)
  }

  if (length(params_to_strip) == 0) {
    return(stanfit)
  }

  sim <- stanfit@sim

  # Build regex patterns for matching flattened names (e.g., "Z[1,1]", "Z_norm[1]", "psi[1]")
  strip_patterns <- paste0("^(", paste(params_to_strip, collapse = "|"), ")(\\[|$)")

  # Filter fnames_oi (flattened parameter names)
  keep_fnames <- !grepl(strip_patterns, sim$fnames_oi)
  sim$fnames_oi <- sim$fnames_oi[keep_fnames]

  # Filter pars_oi and dims_oi (parameter blocks)
  keep_pars <- !(sim$pars_oi %in% params_to_strip)
  sim$pars_oi <- sim$pars_oi[keep_pars]
  sim$dims_oi <- sim$dims_oi[keep_pars]

  # Filter samples from each chain
  for (i in seq_along(sim$samples)) {
    chain_names <- names(sim$samples[[i]])
    keep_samples <- !grepl(strip_patterns, chain_names)
    sim$samples[[i]] <- sim$samples[[i]][keep_samples]
  }

  # Update n_flatnames
  sim$n_flatnames <- length(sim$fnames_oi)

  stanfit@sim <- sim
  return(stanfit)
}
