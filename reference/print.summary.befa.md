# Print BEFA Summary

Displays the summary tables for a fitted BEFA model in a formatted
console output. Includes factor loadings, fit indices, and reliability.

## Usage

``` r
# S3 method for class 'summary.befa'
print(x, digits = 2, cutoff = NULL, sort = NULL, signif_stars = NULL, ...)
```

## Arguments

- x:

  A `summary.befa` object from
  [`summary.befa()`](https://ricardoreysaez.github.io/BayesEFA/reference/summary.befa.md).

- digits:

  Integer. Number of decimal places to display.

- cutoff:

  Numeric. Override cutoff for hiding small loadings.

- sort:

  Logical. Override sorting by primary factor.

- signif_stars:

  Logical. Override significance star display.

- ...:

  Ignored.
