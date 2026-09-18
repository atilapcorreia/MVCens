#' Generate observations from a registered matrix-variate model
#'
#' Validates the common matrix-variate generation contract and dispatches to
#' the generator registered for `model`. The supported models follow the
#' stochastic constructions used in the accompanying articles: matrix-normal
#' row/column covariance for `MVN`, censoring and missingness for `MVNC`,
#' half-normal skewing for `MVSN` and `MVSNC`, a Gamma scale mixture for
#' `MVST`, row-specific half-normal effects for `MVRSN`, and row-specific
#' exponential effects for `MVREN`. All generated observations use the same
#' `p` by `q` shape as `M`; `Sigma` is the `p` by `p` row covariance or scale
#' component and `Psi` is the `q` by `q` column component. The model name is
#' case-insensitive.
#'
#' The following models are available: `MVN`, `MVNC`, `MVSN`, `MVSNC`, `MVST`,
#' `MVRSN`, `MVREN`, `MVNIG`, and `MVVG`. `MVNC`, `MVSNC`, and `MVNIG` can
#' return censored-data structures; `MVVG` uses the additional `rate`
#' argument; `MVST` uses `nu`; and the high-level `MVREN` interface uses the
#' identified convention `lambda_i = 1`.
#'
#' @param model One registered model name. Matching is case-insensitive.
#' @param n Positive integer giving the number of matrix observations.
#' @param M Numeric `p` by `q` location matrix.
#' @param A Optional `p` by `q` skewness or latent-effect matrix. It is
#'   required by skewed and latent-effect models and ignored by `MVN`.
#' @param Sigma Positive-definite `p` by `p` row covariance or scale matrix.
#' @param Psi Positive-definite `q` by `q` column covariance or scale matrix.
#' @param ... Model-specific generation parameters, such as `cens`, `Ind`,
#'   `nu`, `lambda`, `gamma_tilde`, `rate`, or `return_latent`. For censored
#'   generators, `cc = 1` marks censored entries and `LS` stores their upper
#'   censoring limits. For `MVREN`, the supplied `lambda` must be omitted or a
#'   vector of ones.
#' @return For complete-data models, a numeric array with dimensions
#'   `p` by `q` by `n`. Censored generators return a list containing the
#'   observed array and censoring metadata (`cc` and `LS`), and may include
#'   the uncensored array for audit purposes.
#' @details Random-number generation is delegated to the model specification;
#'   this function does not fit parameters or alter the supplied scientific
#'   model. Set the R random seed before calling it when reproducibility is
#'   required. The returned object should be treated as generated data, not as
#'   evidence that a model has been fitted or validated for a real data set.
#' @examples
#' set.seed(123)
#' M <- matrix(0, 3, 4)
#' x <- mv_random(
#'   "MVN", n = 5, M = M,
#'   Sigma = diag(3), Psi = diag(4)
#' )
#' dim(x)
#' @family MVCens main interface
#' @seealso [mv_fit()], [mv_monte_carlo()]
#' @export
mv_random <- function(model, n, M, A = NULL, Sigma, Psi, ...) {
  spec <- get_model_spec(model)
  n <- validate_positive_integer(n, "n")
  spec$validate(n = n, M = M, A = A, Sigma = Sigma, Psi = Psi,
                mode = "generate", ...)
  spec$generate(n = n, M = M, A = A, Sigma = Sigma, Psi = Psi, ...)
}
