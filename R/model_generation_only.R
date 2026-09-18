mvcens_spec_mvnig <- function() {
  new_model_spec(
    name = "MVNIG",
    validate = function(M = NULL, A = NULL, Sigma = NULL, Psi = NULL,
                        cens = 0, Ind = 1, gamma_tilde = 2, mode, ...) {
      if (mode == "fit") return(invisible(TRUE))
      validate_model_parameters(M, Sigma, Psi, A, require_A = TRUE)
      if (!is.numeric(cens) || length(cens) != 1L || cens < 0 || cens > 1) {
        stop("'cens' must lie in [0, 1].", call. = FALSE)
      }
      if (!Ind %in% 1:3) stop("'Ind' must be 1, 2, or 3.", call. = FALSE)
      if (!is.numeric(gamma_tilde) || gamma_tilde <= 0) stop("'gamma_tilde' must be positive.")
      invisible(TRUE)
    },
    generate = function(n, M, A, Sigma, Psi, cens = 0, Ind = 1,
                        gamma_tilde = 2, ...) {
      rmvnig(n, cens, Ind, M, Sigma, Psi, A, gamma_tilde)
    },
    fit_available = FALSE
  )
}

mvcens_spec_mvvg <- function() {
  new_model_spec(
    name = "MVVG",
    validate = function(M = NULL, A = NULL, Sigma = NULL, Psi = NULL,
                        rate = 1, mode, ...) {
      if (mode == "fit") return(invisible(TRUE))
      validate_model_parameters(M, Sigma, Psi, A, require_A = TRUE)
      if (!is.numeric(rate) || rate <= 0) stop("'rate' must be positive.")
      invisible(TRUE)
    },
    generate = function(n, M, A, Sigma, Psi, rate = 1, ...) {
      rmvvg(n, M, A, Sigma, Psi, rate)
    },
    fit_available = FALSE
  )
}
