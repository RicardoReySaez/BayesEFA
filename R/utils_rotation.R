#' Align Posterior Samples using RSP with Varimax Rotation
#'
#' Extracts Lambda draws, corrects for sign/permutation indeterminacy using
#' the RSP algorithm with Varimax rotation, and updates the stanfit object
#' with the aligned values.
#'
#' @param stanfit The original stanfit object.
#' @param J Integer. Number of items.
#' @param M Integer. Number of factors.
#' @param rotation Character. Either "varimax" (default) or "none".
#' @param max_iter Integer. Maximum iterations for the RSP algorithm.
#' @param threshold Numeric. Convergence threshold for RSP.
#' @param verbose Logical. Print progress messages.
#'
#' @return A list with stanfit (updated) and rsp_objective.
#' @keywords internal
#' @noRd
update_aligned_samples <- function(stanfit, J, M,
                                   rotation = "varimax",
                                   max_iter = 1000, threshold = 1e-6,
                                   verbose = FALSE) {
  # ─────────────────────────────── #
  #    Extract and Stack Lambda     #
  # ─────────────────────────────── #

  raw_array <- rstan::extract(stanfit, pars = "Lambda", permuted = FALSE, inc_warmup = FALSE)

  dims <- dim(raw_array)
  n_kept <- dims[1]
  n_chains <- dims[2]
  n_params <- dims[3]

  # Stack Chains Sequentially. Dimensions: [Total_Draws x Params]
  S <- n_kept * n_chains
  lambda_stack <- matrix(NA, nrow = S, ncol = n_params)
  for (i in 1:n_chains) {
    row_start <- (i - 1) * n_kept + 1
    row_end <- i * n_kept
    lambda_stack[row_start:row_end, ] <- raw_array[, i, ]
  }

  # ─────────────────────────────────────────────────── #
  #    Apply RSP Alignment with Varimax (if requested)  #
  # ─────────────────────────────────────────────────── #

  if (rotation == "none") {
    # No alignment - return draws as-is
    aligned_matrix <- lambda_stack
    rsp_objective <- NA
  } else {
    # Apply RSP with Varimax rotation
    rsp_out <- rsp_align(
      lambda_draws = lambda_stack,
      n_items = J,
      n_factors = M,
      n_chains = n_chains,
      format = "column_major",
      max_iter = max_iter,
      threshold = threshold,
      add_names = FALSE
    )

    aligned_matrix <- rsp_out$Lambda_hat_mcmc
    rsp_objective <- rsp_out$objective
  }

  # ──────────────────────────────── #
  #    Inject Lambda into stanfit    #
  # ──────────────────────────────── #

  sim <- stanfit@sim

  # Identify parameter names that match "Lambda"
  all_par_names <- names(sim$samples[[1]])
  param_indices <- grep("^Lambda(\\[|$)", all_par_names)

  # Safety check
  if (length(param_indices) != ncol(aligned_matrix)) {
    warning("Mismatch in parameter count during RSP injection. Skipping alignment.")
    return(list(stanfit = stanfit, rsp_objective = NULL))
  }

  # Loop through chains to inject Lambda
  start_row_stack <- 1
  for (i in 1:n_chains) {
    end_row_stack <- start_row_stack + n_kept - 1
    chain_chunk <- aligned_matrix[start_row_stack:end_row_stack, , drop = FALSE]

    for (col_idx in seq_along(param_indices)) {
      par_list_idx <- param_indices[col_idx]
      full_vec <- sim$samples[[i]][[par_list_idx]]
      total_len <- length(full_vec)
      insertion_start <- total_len - n_kept + 1

      sim$samples[[i]][[par_list_idx]][insertion_start:total_len] <- chain_chunk[, col_idx]
    }

    start_row_stack <- end_row_stack + 1
  }

  # Update stanfit
  stanfit@sim <- sim

  return(list(
    stanfit = stanfit,
    rsp_objective = rsp_objective
  ))
}
