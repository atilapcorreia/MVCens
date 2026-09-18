#' MVCens: Matrix-Variate Models with Censoring and Asymmetry
#'
#' MVCens provides a common interface for generating, fitting, and evaluating
#' matrix-variate normal, censored, skewed, and heavy-tailed models. The main
#' user workflow is built around [mv_random()], [mv_fit()], and
#' [mv_monte_carlo()]. Lower-level density, log-likelihood, moment, and
#' reference-parameter functions are also available for direct calculations.
#'
#' @section Main workflow:
#'
#' - Use [mv_random()] to generate matrix-valued observations.
#' - Use [mv_fit()] to estimate a registered model by ECM.
#' - Use [mv_monte_carlo()] to assess an estimator in repeated samples.
#'
#' @section Supported models:
#' The high-level API recognizes `MVN`, `MVNC`, `MVSN`, `MVSNC`, `MVST`,
#' `MVRSN`, and `MVREN`. Generation is additionally available for `MVNIG` and
#' `MVVG`. Censored models use `cc` as the censoring indicator and `LS` as the
#' array of censoring limits.
#'
#' @section Numerical convention:
#' Matrix observations are stored in arrays of dimensions `p` by `q` by `n`.
#' Covariance structures use a `p` by `p` row component (`Sigma`) and a `q` by
#' `q` column component (`Psi`). Consult each function's help page for its
#' model-specific arguments and return value.
#'
#' @seealso [mv_random()], [mv_fit()], [mv_monte_carlo()]
#' @aliases MVCens
#' @docType package
#' @name MVCens-package
"_PACKAGE"
