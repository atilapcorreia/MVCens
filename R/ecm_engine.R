# Model contracts and central execution engine.

initialize_ecm_state <- function(X, max_iter) {
  dimensions <- dim(X)
  p <- dimensions[1L]
  q <- dimensions[2L]

  list(
    p = p,
    q = q,
    n = dimensions[3L],
    loglik = numeric(max_iter),
    criterion = Inf,
    iteration = 0L,
    mu = apply(X, c(1L, 2L), mean),
    A = apply(X, c(1L, 2L), mean),
    Sigma = diag(p),
    Psi = diag(q),
    Vari = kronecker(diag(q), diag(p))
  )
}

ecm_output_diagnostics <- function(loglik_history, criterion,
                                   monotone_tol = 1e-7) {
  loglik_history <- as.numeric(loglik_history)
  differences <- diff(loglik_history)
  monotone_drops <- differences[
    is.finite(differences) & differences < -monotone_tol
  ]

  list(
    loglik = if (length(loglik_history)) utils::tail(loglik_history, 1L) else NA_real_,
    loglik_history = loglik_history,
    criterion = as.numeric(criterion),
    monotone = length(monotone_drops) == 0L,
    monotone_drops = monotone_drops
  )
}

new_model_spec <- function(name, validate, initialize = NULL, e_step = NULL,
                           cm_step = NULL, loglik = NULL, generate = NULL,
                           parameter_count = NULL, fit = NULL,
                           censored = FALSE, fit_available = TRUE) {
  spec <- list(
    name = name,
    validate = validate,
    initialize = initialize,
    e_step = e_step,
    cm_step = cm_step,
    loglik = loglik,
    generate = generate,
    parameter_count = parameter_count,
    fit = fit,
    censored = censored,
    fit_available = fit_available
  )
  class(spec) <- "mvcens_model_spec"
  spec
}

mvcens_model_specs <- function() {
  specs <- list(
    MVN = mvcens_spec_mvn(),
    MVNC = mvcens_spec_mvnc(),
    MVSN = mvcens_spec_mvsn(),
    MVSNC = mvcens_spec_mvsnc(),
    MVST = mvcens_spec_mvst(),
    MVRSN = mvcens_spec_mvrsn(),
    MVREN = mvcens_spec_mvren(),
    MVNIG = mvcens_spec_mvnig(),
    MVVG = mvcens_spec_mvvg()
  )
  specs
}

get_model_spec <- function(model, require_fit = FALSE) {
  specs <- mvcens_model_specs()
  model <- validate_model_name(model, names(specs))
  spec <- specs[[model]]
  if (require_fit && !isTRUE(spec$fit_available)) {
    stop(sprintf("Model '%s' currently supports generation but not fitting.", model), call. = FALSE)
  }
  spec
}

run_ecm_model <- function(spec, X, cc = NULL, LS = NULL,
                          precision = 1e-6, max_iter = 200L, ...) {
  validate_fit_controls(precision, max_iter)
  spec$validate(X = X, cc = cc, LS = LS, mode = "fit", ...)
  if (!is.function(spec$fit)) {
    stop(sprintf("Model '%s' has no fit implementation.", spec$name), call. = FALSE)
  }
  spec$fit(X = X, cc = cc, LS = LS, precision = precision,
           max_iter = max_iter, ...)
}
