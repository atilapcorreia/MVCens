test_that("MVN generation is reproducible and dimensionally correct", {
  truth <- small_truth()
  set.seed(44)
  first <- mv_random("MVN", 4, truth$M, Sigma = truth$Sigma, Psi = truth$Psi)
  set.seed(44)
  second <- mv_random("MVN", 4, truth$M, Sigma = truth$Sigma, Psi = truth$Psi)
  expect_equal(first, second)
  expect_equal(dim(first), c(2, 2, 4))
})

test_that("censored generation returns the complete fitting contract", {
  truth <- small_truth()
  set.seed(45)
  generated <- mv_random(
    "MVNC", 8, truth$M, Sigma = truth$Sigma, Psi = truth$Psi,
    cens = 0.2, Ind = 1
  )
  expect_named(generated, c("X.cens", "cc", "LS"))
  expect_equal(dim(generated$X.cens), c(2, 2, 8))
  expect_identical(dim(generated$cc), dim(generated$X.cens))
  expect_identical(dim(generated$LS), dim(generated$X.cens))
})

test_that("censored generators support the zero and fully censored boundaries", {
  normal <- small_truth()
  skew <- small_truth(skew = TRUE)

  for (case in list(MVNC = normal, MVSNC = skew)) {
    model <- if (is.null(case$A)) "MVNC" else "MVSNC"
    common <- list(
      model = model, n = 8L, M = case$M, A = case$A,
      Sigma = case$Sigma, Psi = case$Psi, Ind = 1
    )

    complete <- do.call(mv_random, c(common, list(cens = 0)))
    expect_true(all(complete$cc == 0L), info = model)
    expect_true(all(complete$LS == 0), info = model)

    interval <- do.call(mv_random, c(common, list(cens = 1)))
    expect_true(all(interval$cc == 1L), info = model)
    expect_true(all(is.finite(interval$X.cens)), info = model)
    expect_true(all(is.finite(interval$LS)), info = model)
    expect_true(all(interval$X.cens <= interval$LS), info = model)

    expect_error(
      do.call(mv_random, utils::modifyList(common, list(cens = 1, Ind = 2))),
      "wholly uninformative",
      info = model
    )
  }
})
