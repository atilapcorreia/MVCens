#' Unified interface for matrix-variate random generation
#'
#' Generates random matrix-valued observations from one of the
#' matrix-variate distributions implemented in the package. The function
#' dispatches the request to the corresponding model-specific random
#' generator and returns its output in a unified format.
#'
#' @param model Character string specifying the model from which observations
#'   are generated. Available options are `"MVSN"`, `"MVST"`, `"MVNIG"`, `"MVVG"`,
#'   `"MVRSN"`, `"MVREN"`, `"MVNC"`, and `"MVSNC"`.
#' @param n Positive integer specifying the number of matrix-valued
#'   observations to be generated.
#' @param M Numeric matrix specifying the location parameter.
#' @param A Optional numeric matrix specifying the skewness, drift, or
#'   latent-effect parameter. It is required for models containing an
#'   asymmetric or mean-shift component, including `"MVSN"`, `"MVST"`,
#'   `"MVNIG"`, `"MVVG"`, `"MVRSN"`, `"MVREN"`, and `"MVSNC"`.
#' @param U Symmetric positive-definite numeric matrix specifying the row
#'   covariance or scale matrix.
#' @param V Symmetric positive-definite numeric matrix specifying the column
#'   covariance or scale matrix.
#' @param lambda Optional positive numeric vector containing the
#'   row-specific exponential rate parameters for `model = "MVREN"`.
#'   Its length must equal `nrow(M)`. When `NULL`, a vector of ones is used.
#' @param return_latent Logical scalar indicating whether latent variables
#'   should be included in the returned object for models whose generators
#'   support this option. The default is `FALSE`.
#' @param cens Optional numeric scalar specifying the desired censoring
#'   proportion for censored models. It must lie between zero and one and is
#'   used only for `model = "MVNC"` and `model = "MVSNC"`.
#' @param Ind Optional integer specifying the censoring or incomplete-data
#'   mechanism used for `model = "MVNC"` and `model = "MVSNC"`.
#'   The admissible values and their interpretation are determined by the
#'   corresponding censored-data generator.
#' @param nu Positive numeric scalar specifying the degrees-of-freedom
#'   parameter for `model = "MVST"`.
#' @param gamma_tilde Positive numeric scalar specifying the shape parameter
#'   of the mixing distribution for `model = "MVNIG"`. The default is `2`.
#' @param rate Positive numeric scalar specifying the exponential rate
#'   parameter for `model = "MVVG"`. The default is `1`.
#'
#' @details
#' The required arguments depend on the selected model:
#'
#' \describe{
#'   \item{`"MVSN"`}{
#'     Requires `n`, `M`, `A`, `U`, and `V`.
#'   }
#'   \item{`"MVST"`}{
#'     Requires `n`, `M`, `A`, `U`, `V`, and `nu`.
#'   }
#'   \item{`"MVNIG"`}{
#'     Requires `n`, `M`, `A`, `U`, and `V`. The argument
#'     `gamma_tilde` controls the corresponding mixing distribution.
#'   }
#'   \item{`"MVVG"`}{
#'     Requires `n`, `M`, `A`, `U`, and `V`. The argument `rate`
#'     controls the corresponding mixing distribution.
#'   }
#'   \item{`"MVRSN"`}{
#'     Requires `n`, `M`, `A`, `U`, and `V`.
#'   }
#'   \item{`"MVREN"`}{
#'     Requires `n`, `M`, `A`, `U`, and `V`. The optional argument
#'     `lambda` specifies the row-specific exponential rates.
#'   }
#'   \item{`"MVNC"`}{
#'     Requires `n`, `M`, `U`, `V`, `cens`, and `Ind`.
#'   }
#'   \item{`"MVSNC"`}{
#'     Requires `n`, `M`, `A`, `U`, `V`, `cens`, and `Ind`.
#'   }
#' }
#'
#' The generated observations are stored in a three-dimensional numeric
#' array whose dimensions are
#' \eqn{p \times q \times n}, where \eqn{p = \mathrm{nrow}(M)} and
#' \eqn{q = \mathrm{ncol}(M)}.
#'
#' For censored models, the function additionally returns the censored
#' observations, censoring indicators, and censoring limits required by
#' [mv_fit()]. The realized censoring proportion may differ slightly from
#' the value supplied through `cens`, depending on the censoring mechanism
#' implemented by the underlying generator.
#'
#' When `return_latent = TRUE`, latent quantities are included only when the
#' selected model-specific generator supports their return. The names and
#' dimensions of these additional components depend on the selected model.
#'
#' @return A model-dependent list returned by the corresponding random
#'   generator. For complete-data models, the returned object contains at
#'   least:
#'
#'   \describe{
#'     \item{`X`}{
#'       A numeric array with dimensions \eqn{p \times q \times n}
#'       containing the generated matrix-valued observations.
#'     }
#'   }
#'
#'   For censored models, the returned object additionally contains:
#'
#'   \describe{
#'     \item{`X.cens`}{
#'       A numeric array containing the censored or incomplete observations.
#'     }
#'     \item{`cc`}{
#'       A censoring-indicator array having the same dimensions as `X`.
#'     }
#'     \item{`LS`}{
#'       An array containing the censoring limits.
#'     }
#'   }
#'
#'   When `return_latent = TRUE` and the selected generator supports latent
#'   output, additional model-specific latent variables are included in the
#'   returned list.
#'
#' @seealso
#' [mv_fit()], [mvren_mean()], [mvren_covariances()],
#' [mvrsn_mean()], [mvrsn_covariances()]
#'
#' @export

mv_random <- function(model,
                      n = NULL,
                      M = NULL,
                      A = NULL,
                      U = NULL,
                      V = NULL,
                      lambda = NULL,
                      return_latent = FALSE,
                      cens = NULL,
                      Ind = NULL,
                      nu = NULL,
                      gamma_tilde = 2,
                      rate = 1) {

  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop("'model' must be a single character string.")
  }

  model <- toupper(model)

  valid_models <- c("MVN", "MVSN", "MVNC", "MVSNC", "MVST", "MVNIG", "MVVG", "MVRSN", "MVREN")

  if (!(model %in% valid_models)) {
    stop(sprintf(
      "'model' must be one of: %s.",
      paste(valid_models, collapse = ", ")
    ))
  }

  if (model == "MVSN") {
    if (is.null(n) || is.null(M) || is.null(A) || is.null(U) || is.null(V)) {
      stop("For model = 'MVSN', you must provide 'n', 'M', 'A', 'U', and 'V'.")
    }

    return(rmvsn(n = n, M = M, A = A, U = U, V = V))
  }

  if (model == "MVNC") {
    if (is.null(n) || is.null(cens) || is.null(M) || is.null(U) || is.null(V)) {
      stop("For model = 'MVNC', you must provide 'n', 'cens', 'M', 'U', and 'V'.")
    }

    return(
      rmatrix_censored(
        n = n,
        cens = cens,
        Ind = Ind,
        M = M,
        U = U,
        V = V,
        A = NULL,
        dist = "Normal"
      )
    )
  }

  if (model == "MVSNC") {
    if (is.null(n) || is.null(cens) || is.null(M) || is.null(A) || is.null(U) || is.null(V)) {
      stop("For model = 'MVSNC', you must provide 'n', 'cens', 'M', 'A', 'U', and 'V'.")
    }

    return(
      rmatrix_censored(
        n = n,
        cens = cens,
        Ind = Ind,
        M = M,
        U = U,
        V = V,
        A = A,
        dist = "SN"
      )
    )
  }

  if (model == "MVST") {
    if (is.null(n) || is.null(M) || is.null(A) || is.null(U) || is.null(V) || is.null(nu)) {
      stop("For model = 'MVST', you must provide 'n', 'M', 'A', 'U', 'V', and 'nu'.")
    }

    return(rmvst(n = n, M = M, A = A, U = U, V = V, nu = nu))
  }

  if (model == "MVNIG") {
    if (is.null(n) || is.null(cens) || is.null(M) || is.null(A) || is.null(U) || is.null(V)) {
      stop("For model = 'MVNIG', you must provide 'n', 'cens', 'M', 'A', 'U', and 'V'.")
    }

    return(
      rmvnig(
        n = n,
        cens = cens,
        Ind = Ind,
        M = M,
        U = U,
        V = V,
        A = A,
        gamma_tilde = gamma_tilde
      )
    )
  }

  if (model == "MVVG") {
    if (is.null(n) || is.null(M) || is.null(A) || is.null(U) || is.null(V)) {
      stop("For model = 'MVVG', you must provide 'n', 'M', 'A', 'U', and 'V'.")
    }

    return(rmvvg(n = n, M = M, A = A, U = U, V = V, rate = rate))
  }

  if (model == "MVRSN") {
    if (is.null(n) || is.null(M) || is.null(A) || is.null(U) || is.null(V)) {
      stop("For model = 'MVRSN', you must provide 'n', 'M', 'A', 'U', and 'V'.")
    }

    return(
      rmvrsn(
        n = n,
        M = M,
        A = A,
        Sigma = U,
        Psi = V,
        seed = NULL,
        return_latent = FALSE
      )
    )
  }

  if (model == "MVREN") {
    if (is.null(n) || is.null(M) || is.null(A) || is.null(U) || is.null(V) || is.null(lambda)) {
      stop(
        paste0(
          "For model = 'MVREN', you must provide ",
          "'n', 'M', 'A', 'U', 'V' and 'lambda'."
        )
      )
    }

    return(
      rmvren(
        n = n,
        M = M,
        A = A,
        Sigma = U,
        Psi = V,
        lambda = lambda,
        return_latent = FALSE
      )
    )
  }

  stop("Unexpected 'model' value.")
}

#' Unified interface for matrix-variate ECM model fitting
#'
#' Fits one of the matrix-variate models implemented in the package using the
#' corresponding non-internal ECM estimation function.
#'
#' @param model Character string specifying the model to be fitted.
#' @param samples Numeric three-dimensional array containing the observed
#' matrix-valued data.
#' @param cc Optional censoring-indicator array for censored models.
#' @param LS Optional array containing censoring limits.
#' @param precision Positive numeric scalar specifying the convergence
#' tolerance.
#' @param MaxIter Positive integer specifying the maximum number of ECM
#' iterations.
#' @param nu Positive numeric scalar specifying the initial degrees of freedom
#' for the MVST model.
#' @param get.nu Logical scalar indicating whether the degrees-of-freedom
#' parameter should be estimated for the MVST model.
#' @param nu_bounds Numeric vector of length two specifying the admissible
#' interval for the degrees-of-freedom parameter in the MVST model.
#' @param lambda_mode Character string specifying how the exponential-rate
#' parameter vector is handled in the MVREN model. Available options are
#' \code{"fixed"}, \code{"unconstrained"}, and
#' \code{"unit_A_rows"}. Under \code{"fixed"}, the rate parameters remain
#' fixed throughout the ECM iterations. Under \code{"unconstrained"}, both
#' the latent-effect matrix and the rate parameters are estimated without
#' imposing a normalization constraint. Under \code{"unit_A_rows"}, each
#' row of the latent-effect matrix is normalized to have unit Euclidean norm
#' after each update, with the corresponding rate parameter rescaled
#' accordingly.
#' @param M_init Optional initial value for the location matrix. When
#' \code{NULL}, an initial value is constructed internally.
#' @param A_init Optional initial value for the latent-effect matrix. When
#' \code{NULL}, an initial value is constructed internally.
#' @param Sigma_init Optional initial value for the row covariance matrix.
#' When \code{NULL}, an initial value is constructed internally.
#' @param Psi_init Optional initial value for the column covariance matrix.
#' When \code{NULL}, an initial value is constructed internally.
#' @param lambda_init Optional positive numeric vector containing the initial
#' exponential-rate parameters for the MVREN model. Its length must equal
#' the number of rows of each matrix-valued observation. When \code{NULL},
#' the initial rate vector is constructed internally.
#' @param q_policy Character string specifying how a numerically non-positive-
#' definite matrix \code{Q} is handled in the MVREN model. Available options
#' are \code{"warn"}, \code{"strict"}, and \code{"regularize"}.
#'
#' @return An object whose class depends on the selected model:
#'   `"MVN.ECM"`, `"MVSN.ECM"`, `"MVST.ECM"`, `"MVRSN.ECM"`,
#'   `"MVREN.ECM"`, `"MVNC.ECM"`, or `"MVSNC.ECM"`.
#'   The returned object is a list containing the fitted model parameters,
#'   convergence information, and model-specific diagnostic quantities.
#'   Depending on the selected model, the returned list contains:
#'
#'  \describe{
#'    \item{M}{Estimated location matrix.}
#'    \item{A}{Estimated skewness or latent-effect matrix, when applicable.}
#'    \item{Sigma}{Estimated row covariance or scale matrix.}
#'    \item{Psi}{Estimated column covariance or scale matrix.}
#'    \item{lambda}{Estimated or fixed exponential-rate vector for MVREN.}
#'    \item{nu}{Estimated or fixed degrees of freedom for MVST.}
#'    \item{loglik}{Final observed-data log-likelihood.}
#'    \item{BIC}{Bayesian information criterion.}
#'    \item{iterations}{Number of ECM iterations performed.}
#'    \item{converged}{Logical value indicating whether the convergence
#'       criterion was satisfied.}
#'    \item{loglik_history}{Observed-data log-likelihood values recorded
#'       throughout the iterations.}
#'   }
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' p <- 3
#' q <- 4
#' n <- 100
#'
#' M <- matrix(0, p, q)
#' A <- matrix(1, p, q)
#' U <- diag(p)
#' V <- diag(q)
#'
#' ## Matrix-variate Normal
#' sim_mvn <- mv_random(
#'   model = "MVN",
#'   n = n,
#'   M = M,
#'   U = U,
#'   V = V
#' )
#'
#' fit_mvn <- mv_fit(
#'   model = "MVN",
#'   samples = sim_mvn$X
#' )
#'
#' ## Matrix-variate Skew-Normal
#' sim_mvsn <- mv_random(
#'   model = "MVSN",
#'   n = n,
#'   M = M,
#'   A = A,
#'   U = U,
#'   V = V
#' )
#'
#' fit_mvsn <- mv_fit(
#'   model = "MVSN",
#'   samples = sim_mvsn$X
#' )
#'
#' ## Matrix-variate Skew-t
#' sim_mvst <- mv_random(
#'   model = "MVST",
#'   n = n,
#'   M = M,
#'   A = A,
#'   U = U,
#'   V = V,
#'   nu = 5
#' )
#'
#' fit_mvst <- mv_fit(
#'   model = "MVST",
#'   samples = sim_mvst$X,
#'   nu = 5
#' )
#'
#' ## Matrix-variate Row Skew-Normal
#' sim_mvrsn <- mv_random(
#'   model = "MVRSN",
#'   n = n,
#'   M = M,
#'   A = A,
#'   U = U,
#'   V = V
#' )
#'
#' fit_mvrsn <- mv_fit(
#'   model = "MVRSN",
#'   samples = sim_mvrsn$X
#' )
#'
#' ## Matrix-variate Censored Normal
#' cens <- 0.15
#'
#' sim_mvnc <- mv_random(
#'   model = "MVNC",
#'   n = n,
#'   cens = cens,
#'   M = M,
#'   U = U,
#'   V = V
#' )
#'
#' fit_mvnc <- mv_fit(
#'   model = "MVNC",
#'   samples = sim_mvnc$X.cens,
#'   cc = sim_mvnc$cc,
#'   LS = sim_mvnc$LS
#' )
#'
#' ## Matrix-variate Censored Skew-Normal
#' sim_mvsnc <- mv_random(
#'   model = "MVSNC",
#'   n = n,
#'   cens = cens,
#'   M = M,
#'   A = A,
#'   U = U,
#'   V = V
#' )
#'
#' fit_mvsnc <- mv_fit(
#'   model = "MVSNC",
#'   samples = sim_mvsnc$X.cens,
#'   cc = sim_mvsnc$cc,
#'   LS = sim_mvsnc$LS
#' )
#' }

mv_fit <- function(model,
                   samples = NULL,
                   cc = NULL,
                   LS = NULL,
                   precision = 1e-6,
                   MaxIter = 200,
                   nu = 4,
                   get.nu = TRUE,
                   nu_bounds = c(0.05, 150),
                   lambda_mode = c("fixed", "unconstrained", "unit_A_rows"),
                   M_init = NULL,
                   A_init = NULL,
                   Sigma_init = NULL,
                   Psi_init = NULL,
                   lambda_init = NULL,
                   q_policy = c("warn", "strict", "regularize")) {

  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop("'model' must be a single character string.")
  }

  model <- toupper(model)
  valid_models <- c("MVN", "MVNC", "MVSN", "MVSNC", "MVST", "MVRSN", "MVREN")

  if (!(model %in% valid_models)) {
    stop(sprintf(
      "'model' must be one of: %s.",
      paste(valid_models, collapse = ", ")
    ))
  }

  if (model == "MVN") {
    if (is.null(samples)) {
      stop("For model = 'MVN', you must provide 'samples'.")
    }

    return(
      mvn_ecm(
        samples = samples,
        precision = precision,
        MaxIter = MaxIter
      )
    )
  }

  if (model == "MVNC") {
    if (is.null(samples) || is.null(cc) || is.null(LS)) {
      stop("For model = 'MVNC', you must provide 'samples', 'cc', and 'LS'.")
    }

    return(
      mvnc_ecm(
        samples = samples,
        cc = cc,
        LS = LS,
        precision = precision,
        MaxIter = MaxIter
      )
    )
  }

  if (model == "MVSN") {
    if (is.null(samples)) {
      stop("For model = 'MVSN', you must provide 'samples'.")
    }

    return(
      mvsn_ecm(
        dados = samples,
        precision = precision,
        MaxIter = MaxIter,
      )
    )
  }

  if (model == "MVSNC") {
    if (is.null(samples) || is.null(cc) || is.null(LS)) {
      stop("For model = 'MVSNC', you must provide 'samples', 'cc', and 'LS'.")
    }

    return(
      mvsnc_ecm(
        dados = samples,
        cc = cc,
        LS = LS,
        precision = precision,
        MaxIter = MaxIter,
      )
    )
  }

  if (model == "MVST") {
    if (is.null(samples)) {
      stop("For model = 'MVST', you must provide 'samples'.")
    }

    return(
      mvst_ecm(
        dados = samples,
        nu = nu,
        precision = precision,
        MaxIter = MaxIter,
        get.nu = get.nu,
        nu_bounds = nu_bounds
      )
    )
  }

  if (model == "MVRSN") {
    if (is.null(samples)) {
      stop("For model = 'MVRSN', you must provide 'samples'.")
    }

    return(
      mvrsn_ecm(
        X = samples,
        max_iter = MaxIter,
        tol = precision,
        M_init = M_init,
        A_init = A_init,
        Sigma_init = Sigma_init,
        Psi_init = Psi_init,
      )
    )
  }

  if (model == "MVREN") {
    if (is.null(samples)) {
      stop(
        "For model = 'MVREN', you must provide 'samples'.",
        call. = FALSE
      )
    }

    lambda_mode <- match.arg(lambda_mode)
    q_policy <- match.arg(q_policy)

    return(
      mvren_ecm(
        X = samples,
        max_iter = MaxIter,
        tol = precision,
        lambda_mode = lambda_mode,
        M_init = M_init,
        A_init = A_init,
        Sigma_init = Sigma_init,
        Psi_init = Psi_init,
        lambda_init = lambda_init,
        q_policy = q_policy
      )
    )
  }

  stop(
    "Unexpected 'model' value.",
    call. = FALSE
  )
}



