# Shared, deterministic Monte Carlo infrastructure.

`%||%` <- function(left, right) if (!is.null(left)) left else right

default_mc_truth <- function(skew = FALSE, nu = NULL) {
  truth <- list(M = matrix(0, 2L, 2L), Sigma = diag(2L), Psi = diag(2L))
  if (skew) truth$A <- matrix(0.2, 2L, 2L)
  if (!is.null(nu)) truth$nu <- nu
  truth
}

# Each task receives its own L'Ecuyer stream before scheduling. Consequently,
# a campaign is reproducible for a fixed seed regardless of worker count.
run_independent_tasks <- function(tasks, worker, seed, workers = NULL) {
  if (!length(tasks)) return(list())
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed)) {
    stop("'seed' must be one integer.", call. = FALSE)
  }
  if (is.null(workers)) workers <- parallel::detectCores(logical = TRUE)
  workers <- max(1L, min(validate_positive_integer(workers, "workers"), length(tasks)))

  old_kind <- RNGkind()
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (old_seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  streams <- vector("list", length(tasks))
  streams[[1L]] <- get(".Random.seed", envir = .GlobalEnv)
  if (length(tasks) > 1L) {
    for (i in 2:length(tasks)) {
      streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
    }
  }
  seeded_tasks <- Map(function(task, stream) list(task = task, stream = stream),
                      tasks, streams)
  execute <- function(item) {
    assign(".Random.seed", item$stream, envir = .GlobalEnv)
    worker(item$task)
  }

  if (workers == 1L || length(tasks) == 1L) return(lapply(seeded_tasks, execute))
  if (.Platform$OS.type == "windows") {
    cluster <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    return(parallel::parLapplyLB(cluster, seeded_tasks, execute))
  }
  parallel::mclapply(seeded_tasks, execute, mc.preschedule = FALSE,
                     mc.cores = workers, mc.set.seed = FALSE)
}

summarize_mc_results <- function(results) {
  groups <- split(results, results$n)
  summary <- do.call(rbind, lapply(groups, function(data) {
    output <- data.frame(
      n = unique(data$n), replications = nrow(data),
      convergence_rate = mean(data$converged, na.rm = TRUE),
      monotone_rate = mean(data$monotone, na.rm = TRUE),
      median_iterations = stats::median(data$iterations, na.rm = TRUE)
    )
    for (column in grep("^err_", names(data), value = TRUE)) {
      values <- data[[column]]
      if (!all(is.na(values))) {
        output[[paste0("median_", column)]] <- stats::median(values, na.rm = TRUE)
      }
    }
    output
  }))
  rownames(summary) <- NULL
  summary
}

mc_result_row <- function(model, task, fit, truth) {
  result <- data.frame(
    model = model, n = task$n, replication = task$replication,
    err_M = NA_real_, err_A = NA_real_, err_Sigma = NA_real_,
    err_Psi = NA_real_, err_nu = NA_real_, err_lambda = NA_real_,
    loglik = NA_real_, BIC = NA_real_, iterations = NA_integer_,
    converged = FALSE, monotone = NA, det_Psi = NA_real_,
    error = NA_character_, stringsAsFactors = FALSE
  )
  if (inherits(fit, "error")) {
    result$error <- conditionMessage(fit)
    return(result)
  }

  result$err_M <- relative_frobenius_error(fit$M %||% fit$mu, truth$M)
  if (!is.null(truth$A) && !is.null(fit$A)) {
    result$err_A <- relative_frobenius_error(fit$A, truth$A)
  }
  result$err_Sigma <- relative_frobenius_error(fit$Sigma, truth$Sigma)
  result$err_Psi <- relative_frobenius_error(fit$Psi, truth$Psi)
  if (!is.null(truth$nu) && !is.null(fit$nu)) result$err_nu <- abs(fit$nu - truth$nu)
  if (!is.null(truth$lambda) && !is.null(fit$lambda)) {
    result$err_lambda <- relative_frobenius_error(fit$lambda, truth$lambda)
  }
  result$loglik <- utils::tail(fit$loglik, 1L)
  result$BIC <- fit$BIC
  result$iterations <- fit$iterations %||% fit$iter
  result$converged <- isTRUE(fit$converged)
  result$monotone <- fit$monotone %||% NA
  result$det_Psi <- det(fit$Psi)
  result
}

run_model_monte_carlo <- function(model, sample_sizes, replications, truth,
                                  precision, max_iter, seed, workers,
                                  cens = NULL, Ind = NULL, verbose = FALSE,
                                  progress_dir = NULL,
                                  ...) {
  spec <- get_model_spec(model, require_fit = TRUE)
  if (!length(sample_sizes)) {
    stop("'sample_sizes' must contain positive integers.",
         call. = FALSE)
  }
  sample_sizes <- vapply(sample_sizes, validate_positive_integer, integer(1L),
                         name = "sample_sizes")
  replications <- validate_positive_integer(replications, "replications")
  validate_fit_controls(precision, max_iter)
  tasks <- unlist(lapply(sample_sizes, function(n) {
    lapply(seq_len(replications), function(replication) {
      list(n = n, replication = replication)
    })
  }), recursive = FALSE)
  extra <- list(...)

  rows <- run_independent_tasks(tasks, function(task) {
    progress_file <- NULL
    progress_callback <- NULL
    if (!is.null(progress_dir) && identical(model, "MVREN")) {
      dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
      progress_file <- file.path(progress_dir, sprintf(
        "active_%d_%d_%d.tsv", task$n, task$replication, Sys.getpid()
      ))
      progress_callback <- function(iteration, max_iter, criterion) {
        line <- paste(
          Sys.getpid(), task$n, task$replication, iteration, max_iter,
          as.numeric(Sys.time()), sep = "\t"
        )
        tmp <- paste0(progress_file, ".tmp")
        writeLines(line, tmp, useBytes = TRUE)
        file.rename(tmp, progress_file)
        invisible(NULL)
      }
      on.exit(unlink(c(progress_file, paste0(progress_file, ".tmp"))), add = TRUE)
    }
    generation_args <- c(list(
      n = task$n, M = truth$M, A = truth$A %||% NULL,
      Sigma = truth$Sigma, Psi = truth$Psi
    ), extra)
    if (!is.null(truth$nu) && is.null(generation_args$nu)) generation_args$nu <- truth$nu
    if (!is.null(truth$lambda) && is.null(generation_args$lambda)) {
      generation_args$lambda <- truth$lambda
    }
    if (isTRUE(spec$censored)) {
      generation_args$cens <- cens
      generation_args$Ind <- Ind
    }

    fit <- tryCatch({
      generated <- do.call(spec$generate, generation_args)
      fit_args <- c(list(
        spec = spec,
        X = if (isTRUE(spec$censored)) generated$X.cens else generated,
        cc = if (isTRUE(spec$censored)) generated$cc else NULL,
        LS = if (isTRUE(spec$censored)) generated$LS else NULL,
        precision = precision, max_iter = max_iter
      ), extra)
      if (!is.null(progress_callback)) {
        fit_args$progress_callback <- progress_callback
      }
      if (model == "MVST" && is.null(fit_args$nu)) fit_args$nu <- truth$nu
      suppressWarnings(do.call(run_ecm_model, fit_args))
    }, error = function(error) error)
    mc_result_row(model, task, fit, truth)
  }, seed = seed, workers = workers)

  results <- do.call(rbind, rows)
  list(results = results, summary = summarize_mc_results(results), truth = truth)
}

#' Run a Monte Carlo campaign for a registered MVCens model
#'
#' Runs repeated generate-and-fit experiments for a registered estimable model.
#' This is intended to reproduce the Monte Carlo style used in the accompanying
#' matrix-variate articles: data are generated from a known matrix-valued truth,
#' the corresponding EM-type estimator is fitted, and estimation error is
#' summarized across replications. Independent replications are distributed
#' across all available logical CPU threads by default. ECM iterations remain
#' sequential, preserving their mathematical update order. Random streams are
#' reproducible across worker counts for a fixed seed.
#'
#' @param model One of `MVN`, `MVNC`, `MVSN`, `MVSNC`, `MVST`, `MVRSN`, or
#'   `MVREN`.
#' @param sample_sizes Positive integer vector of sample sizes.
#' @param replications Positive number of replications per sample size.
#' @param truth Named list containing the data-generating parameters, typically
#'   `M`, `Sigma`, `Psi`, and, for skewed or latent-effect models, `A`.
#' @param workers Number of parallel workers; `NULL` uses all logical CPUs.
#' @param precision ECM convergence tolerance.
#' @param max_iter Maximum ECM iterations.
#' @param seed Integer random seed.
#' @param progress_dir Optional directory used to publish transient per-worker
#'   MVREN ECM iteration progress. Files are removed when each replication
#'   finishes. `NULL` disables progress reporting.
#' @param ... Model-specific options such as `cens` or `Ind`; censored
#'   campaigns use `cc = 1` for censored entries in generated samples.
#' @return A list with per-replication `results`, grouped `summary`, and the
#'   supplied `truth`. The summary is descriptive and does not change model
#'   defaults or fitted estimates.
#' @examples
#' truth <- list(
#'   M = matrix(0, 3, 4),
#'   Sigma = diag(3),
#'   Psi = diag(4)
#' )
#' study <- mv_monte_carlo(
#'   "MVN", sample_sizes = 5, replications = 1, truth = truth,
#'   workers = 1, max_iter = 1, seed = 123
#' )
#' study$summary
#' @family MVCens main interface
#' @export
mv_monte_carlo <- function(model, sample_sizes = c(50, 100, 200, 400),
                           replications = 200, truth, workers = NULL,
                           precision = 1e-6, max_iter = 500L,
                           seed = 123L, progress_dir = NULL, ...) {
  spec <- get_model_spec(model, require_fit = TRUE)
  if (!length(sample_sizes)) {
    stop("'sample_sizes' must contain positive integers.",
         call. = FALSE)
  }
  sample_sizes <- vapply(sample_sizes, validate_positive_integer, integer(1L),
                         name = "sample_sizes")
  replications <- validate_positive_integer(replications, "replications")
  validate_fit_controls(precision, max_iter)
  dots <- list(...)

  do.call(run_model_monte_carlo, c(list(
    model = spec$name, sample_sizes = sample_sizes,
    replications = replications, truth = truth,
    precision = precision, max_iter = max_iter,
    seed = seed, workers = workers,
    cens = dots$cens %||% NULL, Ind = dots$Ind %||% NULL,
    verbose = dots$verbose %||% FALSE,
    progress_dir = progress_dir
  ), dots[setdiff(names(dots), c("cens", "Ind", "verbose"))]))
}
