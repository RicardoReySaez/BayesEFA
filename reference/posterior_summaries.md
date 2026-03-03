# Compute Posterior Summaries

A wrapper around
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
to compute summaries directly from a `befa` object.

## Usage

``` r
posterior_summaries(befa_fit, pars, ...)
```

## Arguments

- befa_fit:

  A `befa` object.

- pars:

  Character vector. Parameter names to summarize.

- ...:

  Additional arguments passed to
  [`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html).

## Value

A data frame of summaries.

## Examples

``` r
if (FALSE) { # \dontrun{
befa_fit <- befa(
  data = HS_data, n_factors = 3, model = "cor",
  iter = 500, chains = 2, seed = 123
)

# Posterior summaries of communalities
posterior_summaries(befa_fit, pars = "h2")
} # }
```
