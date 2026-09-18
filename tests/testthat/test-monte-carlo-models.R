test_that("generic Monte Carlo engine covers every estimable model", {
  base <- small_truth()
  skew <- small_truth(skew = TRUE)
  cases <- list(
    MVN = list(truth = base, extra = list()),
    MVNC = list(truth = base, extra = list(cens = 0.2, Ind = 1)),
    MVSN = list(truth = skew, extra = list()),
    MVSNC = list(truth = skew, extra = list(cens = 0.2, Ind = 1)),
    MVST = list(truth = small_truth(skew = TRUE, nu = 4),
                 extra = list(get.nu = FALSE)),
    MVRSN = list(truth = skew, extra = list()),
    MVREN = list(truth = c(skew, list(lambda = c(1, 1))),
                 extra = list())
  )

  for (model in names(cases)) {
    case <- cases[[model]]
    campaign <- do.call(mv_monte_carlo, c(list(
      model = model, sample_sizes = 8L, replications = 1L,
      truth = case$truth, workers = 1L, max_iter = 1L, seed = 17L
    ), case$extra))
    expect_equal(nrow(campaign$results), 1L, info = model)
    expect_true(is.na(campaign$results$error), info = model)
    expect_true(is.finite(campaign$summary$monotone_rate), info = model)
  }
})

test_that("censored Monte Carlo models handle zero and informative full censoring", {
  cases <- list(MVNC = small_truth(), MVSNC = small_truth(skew = TRUE))
  boundaries <- list(
    zero = list(cens = 0, Ind = 1),
    all_interval = list(cens = 1, Ind = 1),
    all_mixed = list(cens = 1, Ind = 3)
  )

  for (model in names(cases)) {
    for (boundary in names(boundaries)) {
      campaign <- do.call(mv_monte_carlo, c(list(
        model = model, sample_sizes = 8L, replications = 1L,
        truth = cases[[model]], workers = 1L, max_iter = 1L, seed = 501L
      ), boundaries[[boundary]]))

      expect_equal(nrow(campaign$results), 1L, info = paste(model, boundary))
      expect_true(is.na(campaign$results$error), info = paste(model, boundary))
      expect_true(
        all(is.finite(unlist(campaign$results[c("loglik", "BIC", "det_Psi")]))),
        info = paste(model, boundary)
      )
    }
  }
})

test_that("fully missing Monte Carlo samples fail closed inside each task", {
  cases <- list(MVNC = small_truth(), MVSNC = small_truth(skew = TRUE))

  for (model in names(cases)) {
    campaign <- mv_monte_carlo(
      model, sample_sizes = 8L, replications = 2L, truth = cases[[model]],
      workers = 1L, max_iter = 1L, seed = 502L, cens = 1, Ind = 2
    )
    expect_true(all(!is.na(campaign$results$error)), info = model)
    expect_match(
      unique(campaign$results$error), "wholly uninformative",
      all = TRUE, info = model
    )
  }
})
