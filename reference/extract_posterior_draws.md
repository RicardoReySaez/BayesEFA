# Extract Posterior Draws

A wrapper around the `posterior` package to extract draws from a `befa`
object in various formats.

## Usage

``` r
extract_posterior_draws(
  befa_fit,
  pars,
  format = c("matrix", "array", "df", "list", "rvars")
)
```

## Arguments

- befa_fit:

  A `befa` object.

- pars:

  Character vector. Parameter names to filter (passed to
  [`posterior::subset_draws()`](https://mc-stan.org/posterior/reference/subset_draws.html)).

- format:

  Character. Output format: "matrix", "array", "df", "list", or "rvars".

## Value

The posterior draws in the requested format.

## Examples

``` r
if (FALSE) { # \dontrun{
befa_fit <- befa(
  data = HS_data, n_factors = 3, model = "cor",
  iter = 500, chains = 2, seed = 123
)

# Extract Lambda draws as a matrix
lambda_draws <- extract_posterior_draws(befa_fit, pars = "Lambda", format = "matrix")

# Extract as rvars for easy summaries
lambda_rvars <- extract_posterior_draws(befa_fit, pars = "Lambda", format = "rvars")
} # }
```
