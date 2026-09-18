mvcens_spec_mvsn <- function() {
  new_model_spec(
    name = "MVSN",
    validate = function(X = NULL, M = NULL, A = NULL, Sigma = NULL,
                        Psi = NULL, mode, ...) {
      if (mode == "fit") validate_mvn_input(X)
      else validate_model_parameters(M, Sigma, Psi, A, require_A = TRUE)
      invisible(TRUE)
    },
    initialize = function(X, max_iter = 200L, ...) initialize_ecm_state(X, max_iter),
    loglik = loglik_mvsn,
    generate = function(n, M, A, Sigma, Psi, ...) rmvsn(n, M, A, Sigma, Psi),
    parameter_count = function(p, q) model_parameter_count(p, q, skew = TRUE),
    fit = function(X, cc = NULL, LS = NULL, precision, max_iter,
                   epsilon = 1e-8, ...) {
      mvsn_ecm(X, precision = precision, max_iter = max_iter, epsilon = epsilon)
    }
  )
}
