#' Unified interface for matrix-variate random generators
#'
#' @param model Character string specifying the model.
#' @param n Positive integer. Number of observations.
#' @param M Numeric matrix. Location matrix.
#' @param A Numeric matrix. Skewness/drift matrix.
#' @param U Square numeric matrix. Row covariance matrix.
#' @param V Square numeric matrix. Column covariance matrix.
#' @param cens Numeric scalar. Censoring proportion.
#' @param Ind Integer. Type of incomplete-data mechanism.
#' @param nu Positive scalar. Degrees of freedom for the MVST model.
#' @param gamma_tilde Positive scalar. Shape parameter for the MVNIG model.
#' @param rate Positive scalar. Exponential rate for the MVVG model.
#'
#' @return The output returned by the selected random generator.
#'
#' @export
#'
#' @examples
#' M <- matrix(c(1.0, 1.5, 2.0,
#'               0.8, 1.2, 1.7,
#'               1.4, 1.9, 2.3),
#'             nrow = 3, ncol = 3, byrow = TRUE)
#'
#' A <- matrix(c( 0.25, -0.10,  0.35,
#'               -0.20,  0.30,  0.15,
#'                0.40,  0.05, -0.25),
#'             nrow = 3, ncol = 3, byrow = TRUE)
#'
#' U <- matrix(c(1.50, 0.35, 0.20,
#'               0.35, 1.20, 0.30,
#'               0.20, 0.30, 1.10),
#'             nrow = 3, ncol = 3, byrow = TRUE)
#'
#' V <- matrix(c(1.30, 0.25, 0.10,
#'               0.25, 1.10, 0.20,
#'               0.10, 0.20, 1.40),
#'             nrow = 3, ncol = 3, byrow = TRUE)
#'
#' ## Matrix-variate skew-normal
#' x_mvsn <- mv_random(
#'   model = "MVSN",
#'   n = 50,
#'   M = M,
#'   A = A,
#'   U = U,
#'   V = V
#' )
#'
#' ## Censored matrix-variate normal
#' x_mvnc <- mv_random(
#'   model = "MVNC",
#'   n = 50,
#'   cens = 0.15,
#'   Ind = 1,
#'   M = M,
#'   U = U,
#'   V = V
#' )
#'
#' ## Censored matrix-variate skew-normal
#' x_mvsnc <- mv_random(
#'   model = "MVSNC",
#'   n = 50,
#'   cens = 0.15,
#'   Ind = 3,
#'   M = M,
#'   U = U,
#'   V = V,
#'   A = A,
#' )
#'
#' ## Matrix-variate skew-t
#' x_mvst <- mv_random(
#'   model = "MVST",
#'   n = 50,
#'   M = M,
#'   A = A,
#'   U = U,
#'   V = V,
#'   nu = 6
#' )
#'
#' ## Matrix-variate normal-inverse Gaussian
#' x_mvnig <- mv_random(
#'   model = "MVNIG",
#'   n = 50,
#'   cens = 0.10,
#'   Ind = 3,
#'   M = M,
#'   U = U,
#'   V = V,
#'   A = A,
#'   gamma_tilde = 2.5
#' )
#'
#' ## Matrix-variate variance-gamma
#' x_mvvg <- mv_random(
#'   model = "MVVG",
#'   n = 50,
#'   M = M,
#'   A = A,
#'   U = U,
#'   V = V,
#'   rate = 1.25
#' )
#'
#' dim(x_mvsn)
#' names(x_mvnc)
#' names(x_mvsnc)
#' dim(x_mvst)
#' names(x_mvnig)
#' dim(x_mvvg)

mv_random <- function(model,
                      n = NULL,
                      M = NULL,
                      A = NULL,
                      U = NULL,
                      V = NULL,
                      cens = NULL,
                      Ind = 1,
                      nu = NULL,
                      gamma_tilde = 2,
                      rate = 1) {

  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop("'model' must be a single character string.")
  }

  model <- toupper(model)

  valid_models <- c("MVSN", "MVNC", "MVSNC", "MVST", "MVNIG", "MVVG")

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
        A = A,
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

  stop("Unexpected 'model' value.")
}

#' Unified interface for matrix-variate ECM model fitting
#'
#' Fits one of the matrix-variate models implemented in the package using the
#' corresponding non-internal ECM estimation function.
#'
#' @param model Character string specifying the model to be fitted.
#' Supported options are:
#' \itemize{
#'   \item \code{"MVN"}   - Matrix Variate Normal
#'   \item \code{"MVNC"}  - Censored Matrix Variate Normal
#'   \item \code{"MVSN"}  - Matrix Variate Skew-Normal
#'   \item \code{"MVSNC"} - Censored Matrix Variate Skew-Normal
#'   \item \code{"MVST"}  - Matrix Variate Skew-t
#' }
#' @param samples Numeric 3D array containing the observed matrix-valued data.
#' For censored models, this should contain the observed values or lower limits.
#' @param cc Optional 3D indicator array for censored models. Entries equal to
#' `1` indicate censored/missing values and `0` indicate observed values.
#' @param LS Optional 3D array of upper limits for censored models.
#' @param precision Positive scalar. Convergence tolerance.
#' @param MaxIter Positive integer. Maximum number of ECM iterations.
#' @param epsilon Positive scalar used for numerical stabilization in models that
#' require it.
#' @param nu Positive scalar. Initial degrees of freedom for the \code{"MVST"} model.
#' @param get.nu Logical. If `TRUE`, updates \code{nu} during ECM for the
#' \code{"MVST"} model.
#' @param nu_bounds Numeric vector of length two. Search interval for \code{nu}
#' in the \code{"MVST"} model.
#'
#' @return
#' The fitted object returned by the selected ECM routine.
#'
#' @export
#'
#' @examples
#' # out <- mv_fit(
#' #   model = "MVN",
#' #   samples = samples,
#' #   MaxIter = 20
#' # )
mv_fit <- function(model,
                   samples = NULL,
                   cc = NULL,
                   LS = NULL,
                   precision = 1e-7,
                   MaxIter = 100,
                   epsilon = 1e-8,
                   nu = 4,
                   get.nu = TRUE,
                   nu_bounds = c(0.05, 150)) {

  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop("'model' must be a single character string.")
  }

  model <- toupper(model)
  valid_models <- c("MVN", "MVNC", "MVSN", "MVSNC", "MVST")

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
        epsilon = epsilon
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
        epsilon = epsilon
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
        epsilon = epsilon,
        nu_bounds = nu_bounds
      )
    )
  }

  stop("Unexpected 'model' value.")
}

