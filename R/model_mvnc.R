mvcens_spec_mvnc <- function() {
  new_model_spec(
    name = "MVNC",
    validate = function(X = NULL, cc = NULL, LS = NULL, M = NULL,
                        Sigma = NULL, Psi = NULL, cens = NULL, Ind = 1,
                        mode, ...) {
      if (mode == "fit") validate_mvnc_input(X, cc, LS)
      else {
        validate_model_parameters(M, Sigma, Psi)
        if (!is.numeric(cens) || length(cens) != 1L || !is.finite(cens) ||
            cens < 0 || cens > 1) {
          stop("'cens' must lie in [0, 1].", call. = FALSE)
        }
      }
      invisible(TRUE)
    },
    loglik = loglik_mvnc,
    generate = function(n, M, A = NULL, Sigma, Psi, cens, Ind = 1, ...) {
      rmatrix_censored(n, cens, Ind, M, Sigma, Psi, dist = "Normal")
    },
    parameter_count = function(p, q) model_parameter_count(p, q),
    fit = function(X, cc, LS, precision, max_iter, ...) {
      mvnc_ecm(samples = X, cc = cc, LS = LS,
               precision = precision, max_iter = max_iter)
    },
    censored = TRUE
  )
}
