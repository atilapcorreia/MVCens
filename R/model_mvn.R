mvcens_spec_mvn <- function() {
  new_model_spec(
    name = "MVN",
    validate = function(X = NULL, M = NULL, Sigma = NULL, Psi = NULL,
                        mode, ...) {
      if (mode == "fit") validate_mvn_input(X)
      else validate_model_parameters(M, Sigma, Psi)
      invisible(TRUE)
    },
    initialize = function(X, max_iter = 200L, ...) initialize_ecm_state(X, max_iter),
    loglik = loglik_mvn,
    generate = function(n, M, A = NULL, Sigma, Psi, ...) rmvn_matrix(n, M, Sigma, Psi),
    parameter_count = function(p, q) model_parameter_count(p, q),
    fit = function(X, cc = NULL, LS = NULL, precision, max_iter, ...) {
      mvn_ecm(samples = X, precision = precision, max_iter = max_iter)
    }
  )
}
