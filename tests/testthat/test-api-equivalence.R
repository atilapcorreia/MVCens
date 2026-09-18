test_that("new MVN API delegates exactly to the scientific kernel", {
  truth <- small_truth()
  set.seed(46)
  X <- mv_random("MVN", 10, truth$M, Sigma = truth$Sigma, Psi = truth$Psi)
  direct <- suppressWarnings(MVCens:::mvn_ecm(X, precision = 1e-5, max_iter = 3))
  unified <- suppressWarnings(mv_fit("MVN", X, precision = 1e-5, max_iter = 3))
  expect_equal(unified$mu, direct$mu, tolerance = 0)
  expect_equal(unified$Sigma, direct$Sigma, tolerance = 0)
  expect_equal(unified$Psi, direct$Psi, tolerance = 0)
  expect_equal(unified$loglik, direct$loglik, tolerance = 0)
})

test_that("deprecated samples argument maps to X", {
  truth <- small_truth()
  set.seed(47)
  X <- mv_random("MVN", 6, truth$M, Sigma = truth$Sigma, Psi = truth$Psi)
  expect_warning(
    fit <- mv_fit("MVN", samples = X, max_iter = 50),
    "deprecated"
  )
  expect_s3_class(fit, "MVN.ECM")
})

test_that("zero-censoring fits delegate exactly to complete-data models", {
  cases <- list(
    MVNC = list(base = "MVN", truth = small_truth()),
    MVSNC = list(base = "MVSN", truth = small_truth(skew = TRUE))
  )

  for (model in names(cases)) {
    case <- cases[[model]]
    args <- list(
      model = case$base, n = 10L, M = case$truth$M, A = case$truth$A,
      Sigma = case$truth$Sigma, Psi = case$truth$Psi
    )
    set.seed(480 + match(model, names(cases)))
    X <- do.call(mv_random, args)
    cc <- array(0L, dim = dim(X))
    LS <- array(0, dim = dim(X))

    base_fit <- suppressWarnings(mv_fit(
      case$base, X, precision = 1e-6, max_iter = 3L
    ))
    censored_fit <- suppressWarnings(mv_fit(
      model, X, cc = cc, LS = LS, precision = 1e-6, max_iter = 3L
    ))

    expect_equal(censored_fit$mu, base_fit$mu, tolerance = 0, info = model)
    if (!is.null(base_fit$A)) {
      expect_equal(censored_fit$A, base_fit$A, tolerance = 0, info = model)
    }
    expect_equal(censored_fit$Sigma, base_fit$Sigma, tolerance = 0, info = model)
    expect_equal(censored_fit$Psi, base_fit$Psi, tolerance = 0, info = model)
    expect_equal(censored_fit$loglik, base_fit$loglik, tolerance = 0, info = model)
  }
})

test_that("fully censored finite intervals produce finite diagnostics", {
  cases <- list(MVNC = small_truth(), MVSNC = small_truth(skew = TRUE))

  for (model in names(cases)) {
    truth <- cases[[model]]
    set.seed(490 + match(model, names(cases)))
    generated <- mv_random(
      model, 10L, truth$M, A = truth$A,
      Sigma = truth$Sigma, Psi = truth$Psi, cens = 1, Ind = 1
    )
    fit <- suppressWarnings(mv_fit(
      model, generated$X.cens, cc = generated$cc, LS = generated$LS,
      max_iter = 2L
    ))

    expect_true(all(is.finite(c(fit$loglik, fit$BIC))), info = model)
    expect_true(all(is.finite(c(fit$Sigma, fit$Psi))), info = model)
  }
})

test_that("wholly missing samples fail closed before ECM iteration", {
  dimensions <- c(2L, 2L, 5L)
  lower <- array(-Inf, dim = dimensions)
  upper <- array(Inf, dim = dimensions)
  cc <- array(1L, dim = dimensions)

  for (model in c("MVNC", "MVSNC")) {
    expect_error(
      mv_fit(model, lower, cc = cc, LS = upper, max_iter = 2L),
      "completamente ausentes",
      info = model
    )
  }
})

test_that("MVN and MVNC preserve matrix observations with singleton axes", {
  cases <- list(
    `1x1` = list(M = matrix(0, 1, 1), Sigma = matrix(1, 1, 1), Psi = matrix(1, 1, 1)),
    `1x2` = list(M = matrix(c(0.1, 0.2), 1, 2), Sigma = matrix(1, 1, 1), Psi = diag(2)),
    `2x1` = list(M = matrix(c(0.1, 0.2), 2, 1), Sigma = diag(2), Psi = matrix(1, 1, 1))
  )

  for (shape in names(cases)) {
    truth <- cases[[shape]]
    for (model in c("MVN", "MVNC")) {
      set.seed(700L + match(shape, names(cases)))
      extra <- if (model == "MVNC") list(cens = 0.2, Ind = 1) else list()
      generated <- do.call(mv_random, c(list(
        model = model, n = 8L, M = truth$M,
        Sigma = truth$Sigma, Psi = truth$Psi
      ), extra))
      fit <- if (model == "MVN") {
        suppressWarnings(mv_fit(model, generated, max_iter = 2L))
      } else {
        suppressWarnings(mv_fit(
          model, generated$X.cens, cc = generated$cc,
          LS = generated$LS, max_iter = 2L
        ))
      }
      expect_s3_class(fit, paste0(model, ".ECM"))
    }
  }
})

test_that("MVN and MVNC reject non-array sample inputs", {
  truth <- small_truth()
  expect_error(
    mv_fit("MVN", matrix(1:4, 2, 2), max_iter = 1L),
    "array numerico 3D"
  )
  expect_error(
    mv_fit("MVNC", matrix(1:4, 2, 2), cc = matrix(0, 2, 2),
            LS = matrix(1, 2, 2), max_iter = 1L),
    "array numerico 3D"
  )
})
