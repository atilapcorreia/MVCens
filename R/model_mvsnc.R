mvcens_spec_mvsnc <- function() {
  new_model_spec(
    name = "MVSNC",
    validate = function(X = NULL, cc = NULL, LS = NULL, M = NULL, A = NULL,
                        Sigma = NULL, Psi = NULL, cens = NULL, Ind = 1,
                        mode, ...) {
      if (mode == "fit") validate_mvnc_input(X, cc, LS)
      else {
        validate_model_parameters(M, Sigma, Psi, A, require_A = TRUE)
        if (!is.numeric(cens) || length(cens) != 1L || !is.finite(cens) ||
            cens < 0 || cens > 1) {
          stop("'cens' must lie in [0, 1].", call. = FALSE)
        }
      }
      invisible(TRUE)
    },
    loglik = loglik_mvsnc,
    generate = function(n, M, A, Sigma, Psi, cens, Ind = 1, ...) {
      rmatrix_censored(n, cens, Ind, M, Sigma, Psi, A = A, dist = "SN")
    },
    parameter_count = function(p, q) model_parameter_count(p, q, skew = TRUE),
    fit = function(X, cc, LS, precision, max_iter, epsilon = 1e-8, ...) {
      mvsnc_ecm(X, cc, LS, precision = precision,
                max_iter = max_iter, epsilon = epsilon)
    },
    censored = TRUE
  )
}
