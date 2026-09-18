# Shared validation contracts for the modular MVCens API.

validate_model_name <- function(model, available) {
  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop("'model' must be one non-missing character string.", call. = FALSE)
  }
  normalized <- toupper(model)
  if (!normalized %in% available) {
    stop(sprintf("'model' must be one of: %s.", paste(available, collapse = ", ")), call. = FALSE)
  }
  normalized
}

validate_positive_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0 || value != as.integer(value)) {
    stop(sprintf("'%s' must be a positive integer.", name), call. = FALSE)
  }
  as.integer(value)
}

validate_integer_at_least <- function(value, name, minimum) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < minimum || value != floor(value)) {
    stop(sprintf("'%s' must be an integer greater than or equal to %s.",
                 name, minimum), call. = FALSE)
  }
  as.integer(value)
}

validate_fit_controls <- function(precision, max_iter) {
  if (!is.numeric(precision) || length(precision) != 1L ||
      !is.finite(precision) || precision <= 0) {
    stop("'precision' must be a positive finite scalar.", call. = FALSE)
  }
  validate_positive_integer(max_iter, "max_iter")
  invisible(TRUE)
}

validate_model_parameters <- function(M, Sigma, Psi, A = NULL,
                                      require_A = FALSE) {
  validate_matrix_generator_inputs(
    n = 1L, M = M, Sigma = Sigma, Psi = Psi,
    A = A, require_A = require_A
  )
  invisible(TRUE)
}
