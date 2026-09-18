#' Fit a registered matrix-variate model by ECM
#'
#' Validates the common fit contract and dispatches to the model-specific
#' EM-type implementation described by the corresponding model article. `X`
#' must contain `n` matrix observations stored as a `p` by `q` by `n` array.
#' Matrix-normal and skew-normal censored fits use ECM updates built around
#' truncated moments, `MVST` uses an ECME-style heavy-tailed skew model, and
#' `MVRSN`/`MVREN` use row-specific latent-variable ECM kernels. The model name
#' is case-insensitive and the scientific update equations remain in the model
#' kernel.
#'
#' Models currently fitted by this interface are `MVN`, `MVNC`, `MVSN`,
#' `MVSNC`, `MVST`, `MVRSN`, and `MVREN`. `MVNIG` and `MVVG` are generation-only
#' models in this package and fail closed with an informative error when passed
#' to `mv_fit()`.
#'
#' @param model One estimable model name. Matching is case-insensitive.
#' @param X Numeric array of dimensions `p` by `q` by `n`, containing the
#'   matrix-valued observations. The deprecated `samples` name may be supplied
#'   through `...` for backward compatibility.
#' @param cc Optional censoring indicator array with the same dimensions as
#'   `X`; `1` marks censored/missing entries and `0` marks observed entries.
#' @param LS Optional array of upper censoring limits, with the same dimensions
#'   as `X`. Required by censored models together with `cc`.
#' @param precision Positive finite ECM convergence tolerance.
#' @param max_iter Positive integer maximum number of ECM iterations.
#' @param ... Model-specific controls and initial values, for example `nu`,
#'   `M_init`, `A_init`, `Sigma_init`, `Psi_init`, `normalize_Psi`, or
#'   `verbose`.
#' @return A model-specific fit object containing estimated parameters and,
#'   where implemented, the likelihood path, iteration count, convergence
#'   flag, BIC, and numerical diagnostics.
#' @details ECM iterations are kept sequential because their update order is
#'   part of each model's scientific specification. The function is
#'   observational with respect to the input data: it returns a fit object and
#'   does not modify `X`, `cc`, or `LS` in place. Generation-only models fail
#'   closed instead of silently switching to another estimator.
#' @examples
#' set.seed(123)
#' M <- matrix(0, 3, 4)
#' x <- mv_random(
#'   "MVN", n = 8, M = M,
#'   Sigma = diag(3), Psi = diag(4)
#' )
#' fit <- mv_fit("MVN", X = x, max_iter = 2)
#' names(fit)
#' @family MVCens main interface
#' @seealso [mv_random()], [mv_monte_carlo()]
#' @export
mv_fit <- function(model, X = NULL, cc = NULL, LS = NULL,
                   precision = 1e-6, max_iter = 500L, ...) {
  dots <- list(...)
  if (is.null(X) && !is.null(dots$samples)) {
    X <- dots$samples
    dots$samples <- NULL
    warning("'samples' is deprecated; use 'X'.", call. = FALSE)
  }
  if (is.null(X)) stop("'X' must be supplied.", call. = FALSE)
  spec <- get_model_spec(model, require_fit = TRUE)
  do.call(run_ecm_model, c(list(
    spec = spec, X = X, cc = cc, LS = LS,
    precision = precision, max_iter = max_iter
  ), dots))
}
