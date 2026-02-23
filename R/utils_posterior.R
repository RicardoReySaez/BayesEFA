#' Extract Posterior Draws
#'
#' A wrapper around the `posterior` package to extract draws from a `befa`
#' object in various formats.
#'
#' @param befa_fit A `befa` object.
#' @param pars Character vector. Parameter names to filter (passed to [posterior::subset_draws()]).
#' @param format Character. Output format: "matrix", "array", "df", "list", or "rvars".
#'
#' @return The posterior draws in the requested format.
#'
#' @examples
#' \dontrun{
#' befa_fit <- befa(
#'   data = HS_data, n_factors = 3, model = "cor",
#'   iter = 500, chains = 2, seed = 123
#' )
#'
#' # Extract Lambda draws as a matrix
#' lambda_draws <- extract_posterior_draws(befa_fit, pars = "Lambda", format = "matrix")
#'
#' # Extract as rvars for easy summaries
#' lambda_rvars <- extract_posterior_draws(befa_fit, pars = "Lambda", format = "rvars")
#' }
#'
#' @importFrom posterior as_draws_matrix as_draws_array as_draws_df as_draws_list as_draws_rvars subset_draws
#' @export
extract_posterior_draws <- function(befa_fit, pars, format = c("matrix", "array", "df", "list", "rvars")) {
  valid_classes <- c("befa")
  if (!inherits(befa_fit, valid_classes)) {
    stop("Input must be of class 'befa'.", call. = FALSE)
  }

  format <- match.arg(format)

  draws <- switch(format,
    "matrix" = posterior::as_draws_matrix(befa_fit$stanfit),
    "array"  = posterior::as_draws_array(befa_fit$stanfit),
    "df"     = posterior::as_draws_df(befa_fit$stanfit),
    "list"   = posterior::as_draws_list(befa_fit$stanfit),
    "rvars"  = posterior::as_draws_rvars(befa_fit$stanfit)
  )

  posterior::subset_draws(draws, variable = pars)
}

#' Compute Posterior Summaries
#'
#' A wrapper around [posterior::summarise_draws()] to compute summaries directly
#' from a `befa` object.
#'
#' @param befa_fit A `befa` object.
#' @param pars Character vector. Parameter names to summarize.
#' @param ... Additional arguments passed to [posterior::summarise_draws()].
#'
#' @return A data frame of summaries.
#'
#' @examples
#' \dontrun{
#' befa_fit <- befa(
#'   data = HS_data, n_factors = 3, model = "cor",
#'   iter = 500, chains = 2, seed = 123
#' )
#'
#' # Posterior summaries of communalities
#' posterior_summaries(befa_fit, pars = "h2")
#' }
#'
#' @importFrom posterior as_draws summarise_draws subset_draws
#' @export
posterior_summaries <- function(befa_fit, pars, ...) {
  valid_classes <- c("befa")
  if (!inherits(befa_fit, valid_classes)) {
    stop("Input must be of class 'befa'.", call. = FALSE)
  }

  posterior::summarise_draws(
    posterior::subset_draws(
      posterior::as_draws(befa_fit$stanfit),
      variable = pars
    ),
    ...
  )
}
