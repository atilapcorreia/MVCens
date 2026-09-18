test_that("ECM fits expose the pertinent MVREN-style output contract", {
  base <- small_truth()
  skew <- small_truth(skew = TRUE)
  cases <- list(
    MVN = list(truth = base, generation = list(), fit = list()),
    MVNC = list(truth = base, generation = list(cens = 0.2, Ind = 1), fit = list()),
    MVSN = list(truth = skew, generation = list(), fit = list()),
    MVSNC = list(truth = skew, generation = list(cens = 0.2, Ind = 1), fit = list()),
    MVST = list(
      truth = small_truth(skew = TRUE, nu = 4),
      generation = list(nu = 4), fit = list(nu = 4, get.nu = FALSE)
    )
  )
  common <- c(
    "M", "mu", "Sigma", "Psi", "loglik", "loglik_history", "BIC",
    "iterations", "iter", "converged", "criterion", "monotone",
    "monotone_drops", "normalize_Psi", "npar"
  )

  for (model in names(cases)) {
    case <- cases[[model]]
    set.seed(600L + match(model, names(cases)))
    generated <- do.call(mv_random, c(list(
      model = model, n = 10L, M = case$truth$M, A = case$truth$A,
      Sigma = case$truth$Sigma, Psi = case$truth$Psi
    ), case$generation))

    fit_args <- list(
      model = model,
      X = if (model %in% c("MVNC", "MVSNC")) generated$X.cens else generated,
      cc = if (model %in% c("MVNC", "MVSNC")) generated$cc else NULL,
      LS = if (model %in% c("MVNC", "MVSNC")) generated$LS else NULL,
      max_iter = 3L
    )
    fit <- suppressWarnings(do.call(mv_fit, c(fit_args, case$fit)))

    expect_true(all(common %in% names(fit)), info = model)
    expect_identical(fit$mu, fit$M, info = model)
    expect_identical(fit$iterations, fit$iter, info = model)
    expect_length(fit$loglik, 1L)
    expect_length(fit$loglik_history, fit$iterations)
    expect_identical(fit$loglik, utils::tail(fit$loglik_history, 1L), info = model)
    expect_true(is.numeric(fit$criterion) && length(fit$criterion) == 1L, info = model)
    expect_true(is.logical(fit$monotone) && length(fit$monotone) == 1L, info = model)
    expect_equal(
      fit$monotone_drops,
      diff(fit$loglik_history)[diff(fit$loglik_history) < -1e-7],
      tolerance = 0, info = model
    )
    expect_identical(fit$normalize_Psi, TRUE, info = model)
    expect_true(is.numeric(fit$npar) && length(fit$npar) == 1L, info = model)

    if (model %in% c("MVSN", "MVSNC", "MVST")) {
      expect_true("A" %in% names(fit), info = model)
    }
    if (model == "MVNC") expect_true("dadosPred" %in% names(fit))
    if (model == "MVST") expect_true("nu" %in% names(fit))
  }
})

test_that("legacy output aliases remain available", {
  truth <- small_truth()
  set.seed(620)
  X <- mv_random("MVN", 10L, truth$M, Sigma = truth$Sigma, Psi = truth$Psi)
  fit <- suppressWarnings(mv_fit("MVN", X, max_iter = 2L))

  expect_identical(fit$mu, fit$M)
  expect_identical(fit$iter, fit$iterations)
})

test_that("MVST fitting rejects degrees of freedom at or below two", {
  truth <- small_truth(skew = TRUE, nu = 4)
  X <- mv_random("MVST", 10L, truth$M, A = truth$A,
                 Sigma = truth$Sigma, Psi = truth$Psi, nu = 4)
  expect_error(
    mv_fit("MVST", X, nu = 2, get.nu = FALSE, max_iter = 2L),
    "greater than 2"
  )
})

test_that("MVST accepts an initial nu of four", {
  truth <- small_truth(skew = TRUE, nu = 4)
  X <- mv_random("MVST", 20L, truth$M, A = truth$A,
                 Sigma = truth$Sigma, Psi = truth$Psi, nu = 4)

  fixed <- mv_fit("MVST", X, nu = 4, get.nu = FALSE, max_iter = 3L)
  expect_identical(fixed$nu, 4)
  expect_true(is.finite(fixed$loglik))

  estimated <- mv_fit(
    "MVST", X, nu = 4, get.nu = TRUE,
    nu_bounds = c(2.01, 150), max_iter = 3L
  )
  expect_true(is.finite(estimated$nu))
  expect_gt(estimated$nu, 2)
})
