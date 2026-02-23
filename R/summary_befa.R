#' Summary and Print Methods for BEFA Objects
#'
#' @description
#' Methods to summarize and display results from Bayesian Exploratory Factor Analysis.
#'
#' @importFrom posterior as_draws_matrix as_draws_array summarise_draws default_convergence_measures quantile2
#' @name summary.befa
NULL

# ──────────────────────────────────────── #
#   Helper: Build Fit Summary Table        #
# ──────────────────────────────────────── #

#' Build the Fit Summary Table from a befa_fitmeasures Object
#'
#' Creates the formatted summary table for fit measures, including
#' LOO rows and Chi2_ppp. Used by both summary.befa and print.befa_fitmeasures.
#'
#' @param fit_obj A `befa_fitmeasures` object.
#' @param probs Quantiles for credible intervals (default: c(0.025, 0.975)).
#' @return A numeric matrix with rows for each fit index and columns
#'   Estimate, SD, Lower_95, Upper_95.
#' @keywords internal
#' @noRd
build_fit_summary_table <- function(fit_obj, probs = c(0.025, 0.975)) {
  # Summarise the draws_array
  bayes_summary <- posterior::summarise_draws(
    fit_obj$fit_draws,
    Estimate = mean,
    SD = stats::sd,
    CI_Low = ~ quantile(.x, probs = probs[1], na.rm = TRUE, names = FALSE),
    CI_High = ~ quantile(.x, probs = probs[2], na.rm = TRUE, names = FALSE)
  )

  # Build matrix
  bayes_tab <- as.matrix(bayes_summary[, -1])

  rownames(bayes_tab) <- bayes_summary$variable

  # Add Chi2_ppp (scalar, not a distribution)
  chi2_ppp_row <- c(fit_obj$details$chi2_ppp, NA, NA, NA)
  bayes_tab <- rbind(
    bayes_tab[1, , drop = FALSE], # Chi2
    Chi2_ppp = chi2_ppp_row,
    bayes_tab[-1, , drop = FALSE] # Rest
  )

  # ELPD summary table (Approximate Normal CI from loo object)
  loo_res <- fit_obj$loo_object
  est_elpd <- loo_res$estimates["elpd_loo", "Estimate"]
  se_elpd <- loo_res$estimates["elpd_loo", "SE"]
  est_looic <- loo_res$estimates["looic", "Estimate"]
  se_looic <- loo_res$estimates["looic", "SE"]
  est_p_loo <- loo_res$estimates["p_loo", "Estimate"]
  se_p_loo <- loo_res$estimates["p_loo", "SE"]

  loo_rows <- rbind(
    ELPD  = c(est_elpd, se_elpd, est_elpd + qnorm(0.025) * se_elpd, est_elpd + qnorm(0.975) * se_elpd),
    LOOIC = c(est_looic, se_looic, est_looic + qnorm(0.025) * se_looic, est_looic + qnorm(0.975) * se_looic),
    p_loo = c(est_p_loo, se_p_loo, est_p_loo + qnorm(0.025) * se_p_loo, est_p_loo + qnorm(0.975) * se_p_loo)
  )

  # Combine
  final_tab <- rbind(bayes_tab, loo_rows)
  colnames(final_tab) <- c("Estimate", "SD", "CI_Low", "CI_High")
  return(final_tab)
}

# ──────────────────────────────────────── #
#   Helper: Build Reliability Table        #
# ──────────────────────────────────────── #

#' Build the Reliability Summary Table from a befa_reliability Object
#'
#' Creates the formatted summary table for reliability, splitting
#' into total and subscale omega coefficients.
#'
#' @param rel_obj A `befa_reliability` object.
#' @param probs Quantiles for credible intervals (default: c(0.025, 0.975)).
#' @return A list with `total` (data frame) and `subscales` (data frame).
#' @keywords internal
#' @noRd
build_reliability_table <- function(rel_obj, probs = c(0.025, 0.975)) {
  omega_summary <- posterior::summarise_draws(
    rel_obj$omega_draws,
    Estimate = mean,
    SD = stats::sd,
    CI_Low = ~ quantile(.x, probs = probs[1], names = FALSE),
    CI_High = ~ quantile(.x, probs = probs[2], names = FALSE)
  )

  total <- as.data.frame(omega_summary[1, -1])
  subscales <- as.data.frame(omega_summary[-1, -1])
  M <- nrow(subscales)
  rownames(subscales) <- paste0("F", 1:M)

  list(total = total, subscales = subscales)
}

# ──────────────────────────────────────── #
#   Helper: Format numeric values          #
# ──────────────────────────────────────── #

#' Format Numeric Values for Display
#' @param n Numeric vector to format.
#' @param digits Number of decimal places.
#' @return Character vector of formatted values.
#' @keywords internal
#' @noRd
.format_numeric <- function(n, digits = 2) {
  sapply(n, function(val) {
    if (is.na(val) || val == "") {
      return("")
    }
    if (is.character(val)) {
      return(val)
    }
    sprintf(paste0("%.", digits, "f"), as.numeric(val))
  })
}

# ──────────────────────────────────────── #
#   Helper: Build Loading Matrix           #
# ──────────────────────────────────────── #

#' Extract and Format Loading Matrix from Posterior Draws
#' @param draws_lambda Matrix of posterior draws for Lambda.
#' @param J Number of items.
#' @param M Number of factors.
#' @param probs Quantiles for credible intervals.
#' @return List with loading matrix, significance matrix, and diagnostics.
#' @keywords internal
#' @noRd
.extract_loading_summary <- function(draws_lambda, J, M, probs) {
  # Point estimates
  lambda_means <- colMeans(draws_lambda)
  lambda_mat <- matrix(lambda_means, nrow = J, ncol = M, byrow = FALSE)

  # Significance (CI excludes 0)
  q_low <- apply(draws_lambda, 2, quantile, probs = probs[1])
  q_high <- apply(draws_lambda, 2, quantile, probs = probs[2])
  sig_vec <- (q_low > 0) | (q_high < 0)
  sig_mat <- matrix(sig_vec, nrow = J, ncol = M, byrow = FALSE)

  # Convergence diagnostics
  summ_stats <- posterior::summarise_draws(draws_lambda, posterior::default_convergence_measures())
  rhat_mat <- matrix(summ_stats$rhat, nrow = J, ncol = M, byrow = FALSE)
  bulk_mat <- matrix(summ_stats$ess_bulk, nrow = J, ncol = M, byrow = FALSE)
  tail_mat <- matrix(summ_stats$ess_tail, nrow = J, ncol = M, byrow = FALSE)

  list(
    lambda = lambda_mat,
    sig_matrix = sig_mat,
    rhat = rhat_mat,
    ess_bulk = bulk_mat,
    ess_tail = tail_mat
  )
}

# ──────────────────────────────────────── #
#   Helper: Format Table 1 (Loadings)      #
# ──────────────────────────────────────── #

#' Format Loading Table for Printing
#' @keywords internal
#' @noRd
.format_loading_table <- function(df, sig_mat, f_cols, cutoff, sort, signif_stars, digits) {
  # Sort by primary factor
  if (sort) {
    primary_factor <- apply(abs(df[, f_cols, drop = FALSE]), 1, which.max)
    max_loading <- apply(abs(df[, f_cols, drop = FALSE]), 1, max)
    ord <- order(primary_factor, -max_loading)
    df <- df[ord, ]
    sig_mat <- sig_mat[ord, ]
  }

  # Format each column
  for (col in names(df)) {
    if (col == "Variable") next
    vals <- df[[col]]

    if (col %in% f_cols) {
      idx <- match(col, f_cols)
      formatted <- .format_numeric(vals, digits)
      if (signif_stars) {
        formatted <- ifelse(sig_mat[, idx], paste0(formatted, "*"), paste0(formatted, " "))
      }
      if (cutoff > 0) {
        formatted[abs(as.numeric(vals)) < cutoff] <- ""
      }
      df[[col]] <- formatted
    } else if (col %in% c("EssBulk", "EssTail")) {
      df[[col]] <- sprintf("%.0f", vals)
    } else {
      df[[col]] <- .format_numeric(vals, digits)
    }
  }

  list(df = df, sig_mat = sig_mat)
}

# ──────────────────────────────────────── #
#   Helper: Print Styled Table             #
# ──────────────────────────────────────── #

#' Print a Styled Table with Box Drawing Characters
#' @keywords internal
#' @noRd
.print_styled_table <- function(df, title = NULL) {
  cn <- colnames(df)
  widths <- sapply(seq_along(df), function(i) max(nchar(cn[i]), max(nchar(df[[i]])))) + 2
  total_width <- sum(widths)
  thick_line <- paste0(rep("\u2017", total_width), collapse = "")
  thin_line <- paste0(rep("\u2014", total_width), collapse = "")

  if (!is.null(title)) cat(title, "\n")
  cat(thick_line, "\n")

  for (i in seq_along(df)) {
    align <- ifelse(i == 1, "left", "right")
    cat(format(cn[i], width = widths[i], justify = align))
  }
  cat("\n", thin_line, "\n", sep = "")

  for (i in 1:nrow(df)) {
    for (j in 1:ncol(df)) {
      align <- ifelse(j == 1, "left", "right")
      cat(format(df[i, j], width = widths[j], justify = align))
    }
    cat("\n")
  }
  cat(thick_line, "\n")
  invisible(total_width)
}

# ──────────────────────────────────────────────────── #
#   Main Function: summary.befa                        #
# ──────────────────────────────────────────────────── #

#' Summarize BEFA Results
#'
#' Computes posterior summaries for a fitted Bayesian EFA model, including
#' point estimates, credible intervals, and convergence diagnostics for factor loadings.
#'
#' @param object A `befa` object returned by [befa()].
#' @param probs Numeric vector of length 2. Quantiles for credible intervals (default: 0.025, 0.975).
#' @param cutoff Numeric. Loadings with absolute value below this threshold are hidden in print output.
#' @param sort Logical. If TRUE, items are sorted by their primary factor loading.
#' @param signif_stars Logical. If TRUE, marks loadings whose CI excludes zero with an asterisk.
#' @param ... Ignored.
#'
#' @details
#' The summary extracts posterior draws for loadings (Lambda), computes means and quantiles,
#' and reports MCMC diagnostics (Rhat, ESS). Communalities (h2) and uniquenesses (u2)
#' are derived from the squared loadings.
#'
#' @return A `summary.befa` object containing:
#'   * `estimates`: Data frame with loadings, h2, u2, and diagnostics per item.
#'   * `sig_matrix`: Logical matrix indicating which loadings have CIs excluding zero.
#'   * `probs`: The quantile probabilities used for credible intervals.
#'   * `phi`: Factor correlation matrix (identity for orthogonal rotations).
#'   * `variance`: List with per-factor and total explained variance proportions.
#'   * `fit_indices`: Fit measures from the `befa` object (if available).
#'   * `reliability`: Reliability estimates from the `befa` object (if available).
#'   * `header`: List with model metadata (model_type, n_factors, N, rotation).
#'   * `print_options`: List of display options (cutoff, sort, signif_stars).
#'   * `tables`: Pre-built display tables for loadings, fit measures, and reliability.
#'
#' @examples
#' \dontrun{
#' befa_fit <- befa(
#'   data = HS_data, n_factors = 3, model = "cor",
#'   iter = 500, chains = 2, seed = 123
#' )
#'
#' # Get and print the summary
#' befa_summary <- summary(befa_fit, sort = TRUE, signif_stars = TRUE)
#' print(befa_summary)
#' }
#'
#' @method summary befa
#' @export
summary.befa <- function(object,
                         probs = c(0.025, 0.975),
                         cutoff = 0,
                         sort = FALSE,
                         signif_stars = FALSE,
                         ...) {
  # ─────────────────────────────────── #
  #    Extract dimensions and labels    #
  # ─────────────────────────────────── #

  J <- object$stan_data$J
  M <- object$n_factors
  var_names <- paste0("Item_", 1:J)
  factor_names <- paste0("F", 1:M)

  # ───────────────────────────── #
  #    Extract posterior draws    #
  # ───────────────────────────── #

  draws_lambda <- posterior::as_draws_matrix(object$stanfit, variable = "Lambda")
  draws_lambda <- draws_lambda[, grep("^Lambda", colnames(draws_lambda))]

  loading_info <- .extract_loading_summary(draws_lambda, J, M, probs)
  lambda_mat <- loading_info$lambda

  # ────────────────────────────── #
  #    Build summary data frame    #
  # ────────────────────────────── #

  h2 <- rowSums(lambda_mat^2)
  u2 <- 1 - h2

  df_summary <- data.frame(
    Variable = var_names,
    as.data.frame(lambda_mat),
    h2 = h2,
    u2 = u2,
    Rhat = apply(loading_info$rhat, 1, max, na.rm = TRUE),
    EssBulk = apply(loading_info$ess_bulk, 1, min, na.rm = TRUE),
    EssTail = apply(loading_info$ess_tail, 1, min, na.rm = TRUE)
  )
  colnames(df_summary) <- c("Variable", factor_names, "h2", "u2", "Rhat", "EssBulk", "EssTail")

  # ──────────────────────────────────── #
  #    Factor correlations & variance    #
  # ──────────────────────────────────── #

  # Varimax and none rotations always have orthogonal factors (identity Phi)
  phi_mat <- diag(M)
  colnames(phi_mat) <- rownames(phi_mat) <- factor_names

  ss_loadings <- colSums(lambda_mat^2)
  prop_var <- ss_loadings / J

  # ───────────────────────── #
  #    Build result object    #
  # ───────────────────────── #

  res <- list(
    estimates = df_summary,
    sig_matrix = loading_info$sig_matrix,
    probs = probs,
    phi = phi_mat,
    variance = list(per_factor = prop_var, total_explained = sum(prop_var)),
    fit_indices = object$fit_indices,
    reliability = object$reliability,
    header = list(
      model_type = object$model_type,
      n_factors = object$n_factors,
      N = object$stan_data$N_complete + object$stan_data$N_incomplete,
      rotation = object$rotation
    ),
    print_options = list(
      cutoff = cutoff,
      sort = sort,
      signif_stars = signif_stars
    )
  )

  # ──────────────────────────── #
  #    Prepare display tables    #
  # ──────────────────────────── #

  res$tables$table1_loadings <- df_summary

  if (!is.null(object$fit_indices)) {
    res$tables$table3_fit_measures <- as.data.frame(
      build_fit_summary_table(object$fit_indices, probs)
    )
  }

  if (!is.null(object$reliability)) {
    rel_tables <- build_reliability_table(object$reliability, probs)
    res$tables$table4_reliability <- rel_tables$subscales
    attr(res$tables$table4_reliability, "omega_total") <- rel_tables$total
  }

  class(res) <- "summary.befa"
  return(res)
}

# ───────────────────────────────────── #
#   Print Method: print.summary.befa    #
# ───────────────────────────────────── #

#' Print BEFA Summary
#'
#' Displays the summary tables for a fitted BEFA model in a formatted console output.
#' Includes factor loadings, fit indices, and reliability.
#'
#' @param x A `summary.befa` object from [summary.befa()].
#' @param digits Integer. Number of decimal places to display.
#' @param cutoff Numeric. Override cutoff for hiding small loadings.
#' @param sort Logical. Override sorting by primary factor.
#' @param signif_stars Logical. Override significance star display.
#' @param ... Ignored.
#'
#' @method print summary.befa
#' @export
print.summary.befa <- function(x, digits = 2, cutoff = NULL, sort = NULL, signif_stars = NULL, ...) {
  # ─────────────────────────── #
  #    Recover print options    #
  # ─────────────────────────── #

  if (is.null(cutoff)) cutoff <- x$print_options$cutoff
  if (is.null(sort)) sort <- x$print_options$sort
  if (is.null(signif_stars)) signif_stars <- x$print_options$signif_stars

  tab_count <- 1

  # ────────────────────────────── #
  #    Table 1: Factor Loadings    #
  # ────────────────────────────── #

  df_t1 <- x$tables$table1_loadings
  f_cols <- colnames(df_t1)[2:(1 + x$header$n_factors)]

  formatted <- .format_loading_table(df_t1, x$sig_matrix, f_cols, cutoff, sort, signif_stars, digits)
  df_t1 <- formatted$df

  cat("\n")
  t1_title <- sprintf("Table %d. Factor Loadings (Pattern Matrix)", tab_count)
  t1_width <- .print_styled_table(df_t1, title = t1_title)
  tab_count <- tab_count + 1

  # Footnotes
  pct_total <- x$variance$total_explained * 100
  notes_t1 <- c(
    sprintf("%s rotation applied.", tolower(x$header$rotation)),
    "Diagnostics show worst-case values across factors (max Rhat, min ESS).",
    sprintf("The %d latent factors accounted for %.1f%% of total variance.", x$header$n_factors, pct_total)
  )
  if (signif_stars) {
    notes_t1 <- c(notes_t1, sprintf("(*) %g%% Credible Interval excludes 0.", (x$probs[2] - x$probs[1]) * 100))
  }
  if (cutoff > 0) {
    notes_t1 <- c(notes_t1, sprintf("Loadings with absolute values < %.2f are hidden.", cutoff))
  }
  cat(paste(strwrap(paste("Note:", paste(notes_t1, collapse = " ")), width = t1_width), collapse = "\n"), "\n")

  # ──────────────────────────────────────────── #
  #    Table 2: Factor Correlations (pending)    #
  # ──────────────────────────────────────────── #

  # if (!is.null(x$tables$table2_correlations)) {
  #   cat("\n")
  #   df_t2 <- x$tables$table2_correlations
  #   df_t2[] <- lapply(df_t2, .format_numeric, digits = digits)
  #   df_t2 <- cbind(Factor = rownames(df_t2), df_t2)
  #   t2_title <- sprintf("Table %d. Factor Correlations", tab_count)
  #   .print_styled_table(df_t2, title = t2_title)
  #   tab_count <- tab_count + 1
  # }

  # ─────────────────────────── #
  #    Table 3: Fit Measures    #
  # ─────────────────────────── #

  if (!is.null(x$tables$table3_fit_measures)) {
    cat("\n")
    df_t3 <- x$tables$table3_fit_measures
    df_t3_fmt <- as.data.frame(lapply(df_t3, .format_numeric, digits = digits))
    df_t3_fmt <- cbind(Index = rownames(df_t3), df_t3_fmt)

    if ("Chi2_ppp" %in% df_t3_fmt$Index) {
      df_t3_fmt[which(df_t3_fmt$Index == "Chi2_ppp"), c("SD", "CI_Low", "CI_High")] <- ""
    }

    t3_title <- sprintf("Table %d. Bayesian Fit Measures", tab_count)
    t3_width <- .print_styled_table(df_t3_fmt, title = t3_title)
    tab_count <- tab_count + 1

    fit_notes <- c("Intervals are 95% Credible Intervals.", "PPP: Posterior Predictive p-value (Ideal > .05).", "p_loo/LOOIC derived from PSIS-LOO.")
    cat(paste(strwrap(paste("Note:", paste(fit_notes, collapse = " ")), width = t3_width), collapse = "\n"), "\n")
  }

  # ────────────────────────── #
  #    Table 4: Reliability    #
  # ────────────────────────── #

  if (!is.null(x$tables$table4_reliability)) {
    cat("\n")
    df_t4 <- x$tables$table4_reliability
    om_total <- attr(df_t4, "omega_total")
    df_t4_print <- as.data.frame(lapply(df_t4, .format_numeric, digits = digits))
    df_t4_print <- cbind(Factor = rownames(df_t4), df_t4_print)

    t4_title <- sprintf("Table %d. Factor Reliability (Coefficient Omega)", tab_count)
    t4_width <- .print_styled_table(df_t4_print, title = t4_title)
    tab_count <- tab_count + 1

    rel_note <- sprintf(
      "Full Scale Omega Total: %s [%s, %s]. Omega coefficients use the full posterior distribution.",
      .format_numeric(om_total$Estimate, digits),
      .format_numeric(om_total$CI_Low, digits),
      .format_numeric(om_total$CI_High, digits)
    )
    cat(paste(strwrap(rel_note, width = t4_width), collapse = "\n"), "\n")
  }

  # ────────────────────────── #
  #    Convergence Warnings    #
  # ────────────────────────── #

  max_rhat <- max(as.numeric(x$tables$table1_loadings$Rhat), na.rm = TRUE)
  if (max_rhat > 1.1) {
    cat(sprintf("\nWarning: Potential convergence issues detected (Max Rhat: %.2f).\n", max_rhat))
  }

  invisible(x)
}
