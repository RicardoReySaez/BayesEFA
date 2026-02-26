#' The 'BayesEFA' package.
#'
#' @description Bayesian Exploratory Factor Analysis using Hamiltonian Monte Carlo
#' via Stan. Posterior draws are post-processed with Varimax rotation and the RSP
#' alignment algorithm to resolve rotational indeterminacy. Provides Bayesian fit
#' indices, reliability estimates (Omega), and factor scores with full uncertainty
#' quantification.
#'
#' @docType package
#' @name BayesEFA-package
#' @aliases BayesEFA
#' @useDynLib BayesEFA, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling extract get_stancode
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom loo loo relative_eff
#' @importFrom posterior as_draws_array
#' @importFrom mvnfast dmvn rmvn
#' @importFrom stats complete.cases dnorm na.omit qnorm quantile rnorm sd
#' @importFrom withr local_seed
#'
#' @references
#' Stan Development Team (2024). RStan: the R interface to Stan. R package. https://mc-stan.org
#'
"_PACKAGE"
