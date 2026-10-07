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
                          precision = 1e-6, max_iter = 200L,
                          epsilon = NULL, nu = NULL, get.nu = NULL,
                          nu_bounds = NULL, normalize_Psi = NULL,
                          M_init = NULL, A_init = NULL,
                          Sigma_init = NULL, Psi_init = NULL,
                          q_policy = NULL, verbose = NULL,
                          eig_floor = NULL, monotone_tol = NULL,
                          progress_callback = NULL) {
  validate_fit_controls(precision, max_iter)

  switch(
    spec$name,
    MVNC = spec$validate(X = X, cc = cc, LS = LS, mode = "fit"),
    MVSNC = spec$validate(X = X, cc = cc, LS = LS, mode = "fit"),
    MVST = spec$validate(X = X, nu = nu, mode = "fit"),
    spec$validate(X = X, mode = "fit")
  )

  if (!is.function(spec$fit)) {
    stop(sprintf("Model '%s' has no fit implementation.", spec$name), call. = FALSE)
  }

  fit_args <- list(
    X = X, cc = cc, LS = LS,
    precision = precision, max_iter = max_iter
  )

  add_if_supplied <- function(name, value) {
    if (!is.null(value)) fit_args[[name]] <<- value
  }

  if (identical(spec$name, "MVSN") || identical(spec$name, "MVSNC")) {
    add_if_supplied("epsilon", epsilon)
  }

  if (identical(spec$name, "MVST")) {
    add_if_supplied("nu", nu)
    add_if_supplied("get.nu", get.nu)
    add_if_supplied("nu_bounds", nu_bounds)
    add_if_supplied("epsilon", epsilon)
  }

  if (identical(spec$name, "MVRSN")) {
    add_if_supplied("normalize_Psi", normalize_Psi)
    add_if_supplied("M_init", M_init)
    add_if_supplied("A_init", A_init)
    add_if_supplied("Sigma_init", Sigma_init)
    add_if_supplied("Psi_init", Psi_init)
    add_if_supplied("verbose", verbose)
    add_if_supplied("eig_floor", eig_floor)
    add_if_supplied("monotone_tol", monotone_tol)
  }

  if (identical(spec$name, "MVREN")) {
    add_if_supplied("M_init", M_init)
    add_if_supplied("A_init", A_init)
    add_if_supplied("Sigma_init", Sigma_init)
    add_if_supplied("Psi_init", Psi_init)
    add_if_supplied("q_policy", q_policy)
    add_if_supplied("verbose", verbose)
    add_if_supplied("progress_callback", progress_callback)
  }

  do.call(spec$fit, fit_args)
}
