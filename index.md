# BayesEFA: Bayesian Exploratory Factor Analysis

The `BayesEFA` package provides a straightforward way to estimate
Bayesian Exploratory Factor Analysis (EFA) models via `Stan`.

While the package supports both `rstan` and `cmdstanr`, the use of
`cmdstanr` is highly recommended as it is the only backend that ensures
compatibility with the latest version of `Stan` and typically offers
faster estimation times. Please install `cmdstanr` by following the
[official guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

## Installation in Windows and Linux

You can install the development version of BayesEFA from
[GitHub](https://github.com/RicardoReySaez/BayesEFA) with:

``` r
# install.packages("devtools")
devtools::install_github("RicardoReySaez/BayesEFA", force = TRUE)
```

Please be a little patient during installation! It will take several
minutes to compile the underlying Stan models. Grab a coffee and don’t
worry if it seems to be taking a while, it’s completely normal.

## Installation in macOS

Since `BayesEFA` includes custom C++ functions, macOS users must
configure the R compilation toolchain. Following the recommendations
from the [latent](https://github.com/Marcosjnez/latent) package, we
suggest using `macrtools` by James Balamuta to automate this setup:

``` r
# 1. Install macrtools
# install.packages("remotes")
remotes::install_github("coatless-mac/macrtools")

# 2. Setup the toolchain and OpenMP support
macrtools::macos_rtools_install()
macrtools::openmp_install()

# 3. Install BayesEFA
devtools::install_github("RicardoReySaez/BayesEFA", force = TRUE)
```

> **Important:** GitHub Actions routinely test the package against the
> latest macOS environments, so keeping your system updated is the best
> way to avoid compilation errors. If errors persist, we highly
> recommend pasting the compilation output into a Large Language Model
> (e.g., ChatGPT, Gemini, Claude) for quick troubleshooting.

## Quick Start

Here is a basic example showing how to fit a Bayesian EFA model:

``` r
library(BayesEFA)

# Fit a 3-factor Bayesian EFA model using the Holzinger & Swineford dataset
fit <- befa(data = HS_data, n_factors = 3)

# Print a comprehensive summary of the model (inspired in psych package output)
summary(befa_fit, cutoff = 0.3, signif_stars = TRUE)

# Table 1. Factor Loadings (Pattern Matrix)
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Variable       F1     F2     F3    h2    u2  Rhat  EssBulk  EssTail
# ———————————————————————————————————————————————————————————————————
# Item_1      0.59*  0.31*         0.47  0.53  1.00     3458     2849
# Item_2      0.46*                0.23  0.77  1.00     4053     2741
# Item_3      0.65*                0.44  0.56  1.00     4146     2274
# Item_4             0.83*         0.71  0.29  1.00     4156     2607
# Item_5             0.86*         0.75  0.25  1.00     4179     2983
# Item_6             0.81*         0.68  0.32  1.00     4291     2817
# Item_7                    0.68*  0.48  0.52  1.00     3343     2081
# Item_8                    0.69*  0.52  0.48  1.00     3614     2593
# Item_9      0.40*         0.49*  0.43  0.57  1.00     3014     2757
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Note: varimax rotation applied. Diagnostics show worst-case values
# across factors (max Rhat, min ESS). The 3 latent factors accounted
# for 52.2% of total variance. (*) 95% Credible Interval excludes 0.
# Loadings with absolute values < 0.30 are hidden.
#
# Table 2. Bayesian Fit Measures
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Index         Estimate     SD    CI_Low   CI_High
# —————————————————————————————————————————————————
# Chi2             47.12   6.98     35.37     62.56
# Chi2_ppp          0.11
# Chi2_Null       918.85   0.00    918.85    918.85
# BRMSEA            0.05   0.02      0.00      0.09
# BGamma            0.99   0.00      0.98      1.00
# Adj_BGamma        0.97   0.02      0.94      1.00
# BMc               0.98   0.01      0.96      1.00
# SRMR              0.05   0.01      0.03      0.06
# BCFI              0.99   0.01      0.97      1.00
# BTLI              0.99   0.02      0.94      1.00
# ELPD          -3416.31  42.49  -3499.60  -3333.02
# LOOIC          6832.62  84.99   6666.05   6999.20
# p_loo            24.97   1.78     21.48     28.46
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Note: Intervals are 95% Credible Intervals. PPP:
# Posterior Predictive p-value (Ideal > .05).
# p_loo/LOOIC derived from PSIS-LOO.
#
# Table 3. Factor Reliability (Coefficient Omega)
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Factor    Estimate    SD  CI_Low  CI_High
# —————————————————————————————————————————
# F1            0.60  0.04    0.52     0.67
# F2            0.73  0.02    0.69     0.76
# F3            0.55  0.04    0.45     0.62
# ‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗
# Full Scale Omega Total: 0.84 [0.82,
# 0.86]. Omega coefficients use the full
# posterior distribution.
```

For a complete guide on how to perform a Bayesian EFA from scratch,
please visit the [Getting
Started](https://ricardoreysaez.github.io/BayesEFA/articles/getting_started.md)
tutorial or explore the comprehensive package website!
