# Public API contract tests for MVCens

# These tests intentionally verify the explicit public signatures of the
# high-level API. The package no longer uses `...` in mv_random(), mv_fit(),
# or mv_monte_carlo(); model-specific arguments are listed explicitly.


test_that("mv_random uses the canonical explicit argument names", {
  expected <- c(
    "model",
    "n",
    "M",
    "A",
    "Sigma",
    "Psi",
    "cens",
    "Ind",
    "nu",
    "lambda",
    "gamma_tilde",
    "rate",
    "return_latent"
  )

  actual <- names(formals(MVCens::mv_random))

  expect_identical(actual, expected)
  expect_false("..." %in% actual)
})


test_that("mv_fit uses the canonical explicit argument names", {
  expected <- c(
    "model",
    "X",
    "cc",
    "LS",
    "precision",
    "max_iter",
    "samples",
    "epsilon",
    "nu",
    "get.nu",
    "nu_bounds",
    "normalize_Psi",
    "M_init",
    "A_init",
    "Sigma_init",
    "Psi_init",
    "q_policy",
    "verbose",
    "eig_floor",
    "monotone_tol",
    "progress_callback"
  )

  actual <- names(formals(MVCens::mv_fit))

  expect_identical(actual, expected)
  expect_false("..." %in% actual)
})


test_that("mv_monte_carlo uses the canonical explicit argument names", {
  expected <- c(
    "model",
    "sample_sizes",
    "replications",
    "truth",
    "workers",
    "precision",
    "max_iter",
    "seed",
    "progress_dir",
    "cens",
    "Ind",
    "verbose",
    "nu",
    "lambda",
    "epsilon",
    "get.nu",
    "nu_bounds",
    "normalize_Psi",
    "M_init",
    "A_init",
    "Sigma_init",
    "Psi_init",
    "q_policy",
    "eig_floor",
    "monotone_tol"
  )

  actual <- names(formals(MVCens::mv_monte_carlo))

  expect_identical(actual, expected)
  expect_false("..." %in% actual)
})
