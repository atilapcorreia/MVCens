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
      output[[paste0("median_", column)]] <- if (all(is.na(values))) {
        NA_real_
      } else {
        stats::median(values, na.rm = TRUE)
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
                                  nu = NULL, lambda = NULL,
                                  epsilon = NULL, get.nu = NULL,
                                  nu_bounds = NULL, normalize_Psi = NULL,
                                  M_init = NULL, A_init = NULL,
                                  Sigma_init = NULL, Psi_init = NULL,
                                  q_policy = NULL, eig_floor = NULL,
                                  monotone_tol = NULL) {
  spec <- get_model_spec(model, require_fit = TRUE)
  if (!length(sample_sizes)) {
    stop("'sample_sizes' must contain positive integers.", call. = FALSE)
  }
  sample_sizes <- vapply(sample_sizes, validate_positive_integer, integer(1L),
                         name = "sample_sizes")
  replications <- validate_positive_integer(replications, "replications")
  validate_fit_controls(precision, max_iter)

  if (is.null(Ind)) Ind <- 1L
  if (isTRUE(spec$censored) && is.null(cens)) {
    stop("'cens' must be supplied for censored Monte Carlo campaigns.",
         call. = FALSE)
  }
  if (identical(model, "MVST") && is.null(nu) && is.null(truth$nu)) {
    stop("'nu' must be supplied either directly or in 'truth' for MVST.",
         call. = FALSE)
  }

  tasks <- unlist(lapply(sample_sizes, function(n) {
    lapply(seq_len(replications), function(replication) {
      list(n = n, replication = replication)
    })
  }), recursive = FALSE)

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

    generation_args <- list(
      n = task$n,
      M = truth$M,
      A = truth$A %||% NULL,
      Sigma = truth$Sigma,
      Psi = truth$Psi
    )

    if (identical(model, "MVNC") || identical(model, "MVSNC")) {
      generation_args$cens <- cens
      generation_args$Ind <- Ind
    }

    if (identical(model, "MVST")) {
      generation_nu <- nu %||% truth$nu
      generation_args$nu <- generation_nu
    }

    if (identical(model, "MVREN")) {
      generation_args$lambda <- lambda %||% truth$lambda
    }

    fit <- tryCatch({
      generated <- do.call(spec$generate, generation_args)
      generated_X <- if (isTRUE(spec$censored)) generated$X.cens else generated
      generated_cc <- if (isTRUE(spec$censored)) generated$cc else NULL
      generated_LS <- if (isTRUE(spec$censored)) generated$LS else NULL

      fit_nu <- if (identical(model, "MVST")) {
        nu %||% truth$nu
      } else {
        NULL
      }

      suppressWarnings(run_ecm_model(
        spec = spec,
        X = generated_X,
        cc = generated_cc,
        LS = generated_LS,
        precision = precision,
        max_iter = max_iter,
        epsilon = epsilon,
        nu = fit_nu,
        get.nu = get.nu,
        nu_bounds = nu_bounds,
        normalize_Psi = normalize_Psi,
        M_init = M_init,
        A_init = A_init,
        Sigma_init = Sigma_init,
        Psi_init = Psi_init,
        q_policy = q_policy,
        verbose = verbose,
        eig_floor = eig_floor,
        monotone_tol = monotone_tol,
        progress_callback = progress_callback
      ))
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
#'   `MVST` may also obtain `nu` from this list and `MVREN` may obtain `lambda`
#'   from it.
#' @param workers Number of parallel workers; `NULL` uses all logical CPUs.
#' @param precision ECM convergence tolerance.
#' @param max_iter Maximum ECM iterations.
#' @param seed Integer random seed.
#' @param progress_dir Optional directory used to publish transient per-worker
#'   MVREN ECM iteration progress. Files are removed when each replication
#'   finishes. `NULL` disables progress reporting.
#' @param cens Censoring proportion for `MVNC` and `MVSNC`. It must be supplied
#'   for censored campaigns.
#' @param Ind Integer censoring/missingness mechanism for `MVNC` and `MVSNC`:
#'   `1` for interval censoring, `2` for missing values, and `3` for a mixture.
#' @param verbose Logical verbosity control passed to models that support it.
#' @param nu Optional degrees of freedom for `MVST`. If `NULL`, the value in
#'   `truth$nu` is used.
#' @param lambda Optional exponential-rate vector used when generating `MVREN`
#'   samples. If `NULL`, `truth$lambda` is used; if that is also `NULL`, the
#'   MVREN generator uses its identified vector of ones.
#' @param epsilon Optional numerical tolerance for `MVSN`, `MVSNC`, and `MVST`.
#' @param get.nu Optional logical control for updating `nu` in `MVST`.
#' @param nu_bounds Optional length-two interval for the MVST `nu` update.
#' @param normalize_Psi Optional logical control for MVRSN normalization of
#'   `Psi`.
#' @param M_init Optional initial location matrix for `MVRSN` and `MVREN`.
#' @param A_init Optional initial latent-effect/skewness matrix for `MVRSN` and
#'   `MVREN`.
#' @param Sigma_init Optional initial row covariance matrix for `MVRSN` and
#'   `MVREN`.
#' @param Psi_init Optional initial column covariance matrix for `MVRSN` and
#'   `MVREN`.
#' @param q_policy Optional MVREN `Q`-matrix handling policy: `"warn"`,
#'   `"strict"`, or `"regularize"`.
#' @param eig_floor Optional positive eigenvalue floor for `MVRSN`.
#' @param monotone_tol Optional MVRSN likelihood-monotonicity tolerance.
#'
#' @return A list with per-replication `results`, grouped `summary`, and the
#'   supplied `truth`. The summary contains one `median_err_*` column for each
#'   error measure present in `results`; when no valid estimate is available
#'   for a measure within a sample-size group, its median is returned as `NA`
#'   rather than omitting the column. The summary is descriptive and does not
#'   change model defaults or fitted estimates.
#'
#' @examples
#' ## Each campaign below uses one replication and one ECM iteration only to
#' ## demonstrate the interface. Increase both values for a real simulation.
#' \donttest{
#' p <- 2
#' q <- 2
#' M <- matrix(0, p, q)
#' A <- matrix(c(0.4, 0.2, 0.1, 0.3), p, q)
#' Sigma <- matrix(c(1.0, 0.2, 0.2, 1.0), p, p)
#' Psi <- matrix(c(1.0, 0.1, 0.1, 1.0), q, q)
#' truth_mvn <- list(M = M, Sigma = Sigma, Psi = Psi)
#' truth_skew <- list(M = M, A = A, Sigma = Sigma, Psi = Psi)
#'
#' ## 1. MVN: one or several sample sizes can be studied.
#' mc_mvn <- mv_monte_carlo(
#'   model = "MVN", sample_sizes = c(20, 30), replications = 1,
#'   truth = truth_mvn, workers = 1,
#'   precision = 1e-6, max_iter = 1, seed = 123
#' )
#' mc_mvn$summary
#'
#' ## 2. MVNC: all three censoring/missingness mechanisms.
#' mc_mvnc_interval <- mv_monte_carlo(
#'   "MVNC", sample_sizes = 20, replications = 1, truth = truth_mvn,
#'   cens = 0.15, Ind = 1, workers = 1, max_iter = 1, seed = 123
#' )
#' mc_mvnc_missing <- mv_monte_carlo(
#'   "MVNC", sample_sizes = 20, replications = 1, truth = truth_mvn,
#'   cens = 0.15, Ind = 2, workers = 1, max_iter = 1, seed = 123
#' )
#' mc_mvnc_mixed <- mv_monte_carlo(
#'   "MVNC", sample_sizes = 20, replications = 1, truth = truth_mvn,
#'   cens = 0.15, Ind = 3, workers = 1, max_iter = 1, seed = 123
#' )
#'
#' ## 3. MVSN: epsilon is forwarded to the skew-normal ECM fit.
#' mc_mvsn <- mv_monte_carlo(
#'   "MVSN", sample_sizes = 20, replications = 1, truth = truth_skew,
#'   epsilon = 1e-8, workers = 1, max_iter = 1, seed = 123
#' )
#'
#' ## 4. MVSNC: all three mechanisms are supported here as well.
#' mc_mvsnc_interval <- mv_monte_carlo(
#'   "MVSNC", sample_sizes = 20, replications = 1, truth = truth_skew,
#'   cens = 0.15, Ind = 1, epsilon = 1e-8,
#'   workers = 1, max_iter = 1, seed = 123
#' )
#' mc_mvsnc_missing <- mv_monte_carlo(
#'   "MVSNC", sample_sizes = 20, replications = 1, truth = truth_skew,
#'   cens = 0.15, Ind = 2, epsilon = 1e-8,
#'   workers = 1, max_iter = 1, seed = 123
#' )
#' mc_mvsnc_mixed <- mv_monte_carlo(
#'   "MVSNC", sample_sizes = 20, replications = 1, truth = truth_skew,
#'   cens = 0.15, Ind = 3, epsilon = 1e-8,
#'   workers = 1, max_iter = 1, seed = 123
#' )
#'
#' ## 5. MVST: nu can come from truth or be supplied directly.
#' truth_mvst <- c(truth_skew, list(nu = 6))
#' mc_mvst_fixed <- mv_monte_carlo(
#'   "MVST", sample_sizes = 20, replications = 1, truth = truth_mvst,
#'   get.nu = FALSE, epsilon = 1e-8,
#'   workers = 1, max_iter = 1, seed = 123
#' )
#' mc_mvst_estimated <- mv_monte_carlo(
#'   "MVST", sample_sizes = 20, replications = 1, truth = truth_mvst,
#'   get.nu = TRUE, nu_bounds = c(2.01, 50), epsilon = 1e-8,
#'   workers = 1, max_iter = 1, seed = 123
#' )
#' mc_mvst_nu_argument <- mv_monte_carlo(
#'   "MVST", sample_sizes = 20, replications = 1, truth = truth_skew,
#'   nu = 6, get.nu = FALSE,
#'   workers = 1, max_iter = 1, seed = 123
#' )
#'
#' ## 6. MVRSN: Psi normalization can be enabled or disabled.
#' mc_mvrsn_normalized <- mv_monte_carlo(
#'   "MVRSN", sample_sizes = 20, replications = 1, truth = truth_skew,
#'   normalize_Psi = TRUE, eig_floor = 1e-8, monotone_tol = 1e-7,
#'   workers = 1, max_iter = 1, seed = 123
#' )
#' mc_mvrsn_custom <- mv_monte_carlo(
#'   "MVRSN", sample_sizes = 20, replications = 1, truth = truth_skew,
#'   normalize_Psi = FALSE,
#'   M_init = M, A_init = A, Sigma_init = Sigma, Psi_init = Psi,
#'   verbose = FALSE, workers = 1, max_iter = 1, seed = 123
#' )
#'
#' ## 7. MVREN: lambda is fixed to one under the identified parameterization.
#' truth_mvren <- c(truth_skew, list(lambda = rep(1, p)))
#' mc_mvren_warn <- mv_monte_carlo(
#'   "MVREN", sample_sizes = 20, replications = 1, truth = truth_mvren,
#'   q_policy = "warn", workers = 1, max_iter = 1, seed = 123
#' )
#' mc_mvren_strict <- mv_monte_carlo(
#'   "MVREN", sample_sizes = 20, replications = 1, truth = truth_mvren,
#'   q_policy = "strict", workers = 1, max_iter = 1, seed = 123
#' )
#' mc_mvren_regularized <- mv_monte_carlo(
#'   "MVREN", sample_sizes = 20, replications = 1, truth = truth_skew,
#'   lambda = rep(1, p), q_policy = "regularize",
#'   M_init = M, A_init = A, Sigma_init = Sigma, Psi_init = Psi,
#'   workers = 1, max_iter = 1, seed = 123
#' )
#'
#' ## Optional MVREN transient progress files.
#' mc_mvren_progress <- mv_monte_carlo(
#'   "MVREN", sample_sizes = 20, replications = 1, truth = truth_mvren,
#'   q_policy = "warn", progress_dir = tempdir(),
#'   workers = 1, max_iter = 1, seed = 123
#' )
#'
#' }
#' @family MVCens main interface
#' @export
mv_monte_carlo <- function(model, sample_sizes = c(50, 100, 200, 400),
                           replications = 200, truth, workers = NULL,
                           precision = 1e-6, max_iter = 500L,
                           seed = 123L, progress_dir = NULL,
                           cens = NULL, Ind = 1L, verbose = FALSE,
                           nu = NULL, lambda = NULL,
                           epsilon = NULL, get.nu = NULL,
                           nu_bounds = NULL, normalize_Psi = NULL,
                           M_init = NULL, A_init = NULL,
                           Sigma_init = NULL, Psi_init = NULL,
                           q_policy = NULL, eig_floor = NULL,
                           monotone_tol = NULL) {
  spec <- get_model_spec(model, require_fit = TRUE)
  if (!length(sample_sizes)) {
    stop("'sample_sizes' must contain positive integers.", call. = FALSE)
  }
  sample_sizes <- vapply(sample_sizes, validate_positive_integer, integer(1L),
                         name = "sample_sizes")
  replications <- validate_positive_integer(replications, "replications")
  validate_fit_controls(precision, max_iter)

  run_model_monte_carlo(
    model = spec$name,
    sample_sizes = sample_sizes,
    replications = replications,
    truth = truth,
    precision = precision,
    max_iter = max_iter,
    seed = seed,
    workers = workers,
    cens = cens,
    Ind = Ind,
    verbose = verbose,
    progress_dir = progress_dir,
    nu = nu,
    lambda = lambda,
    epsilon = epsilon,
    get.nu = get.nu,
    nu_bounds = nu_bounds,
    normalize_Psi = normalize_Psi,
    M_init = M_init,
    A_init = A_init,
    Sigma_init = Sigma_init,
    Psi_init = Psi_init,
    q_policy = q_policy,
    eig_floor = eig_floor,
    monotone_tol = monotone_tol
  )
}
