test_that("public API uses canonical argument names", {
  expect_equal(
    names(formals(mv_random))[1:6],
    c("model", "n", "M", "A", "Sigma", "Psi")
  )
  expect_equal(
    names(formals(mv_fit))[1:7],
    c("model", "X", "cc", "LS", "precision", "max_iter", "...")
  )
  expect_equal(
    names(formals(mv_monte_carlo))[1:5],
    c("model", "sample_sizes", "replications", "truth", "workers")
  )
})

test_that("model-specific Monte Carlo wrappers are internal", {
  wrappers <- c(
    "mvn_monte_carlo", "mvnc_monte_carlo", "mvsn_monte_carlo",
    "mvsnc_monte_carlo", "mvst_monte_carlo", "mvrsn_monte_carlo",
    "mvren_monte_carlo"
  )
  exports <- getNamespaceExports("MVCens")

  expect_true("mv_monte_carlo" %in% exports)
  expect_false(any(wrappers %in% exports))
  expect_true(all(vapply(
    wrappers, exists, logical(1), envir = asNamespace("MVCens"),
    inherits = FALSE
  )))
})

test_that("every exported function has a generated help topic", {
  exports <- getNamespaceExports("MVCens")
  source_man <- test_path("..", "..", "man")
  rd_objects <- if (dir.exists(source_man)) {
    lapply(list.files(source_man, pattern = "[.]Rd$", full.names = TRUE),
           tools::parse_Rd)
  } else {
    unname(tools::Rd_db("MVCens"))
  }
  aliases <- unique(unlist(lapply(rd_objects, function(rd) {
    alias_nodes <- Filter(function(node) {
      identical(attr(node, "Rd_tag"), "\\alias")
    }, rd)
    vapply(alias_nodes, paste, character(1), collapse = "")
  })))

  expect_true(all(exports %in% aliases),
              info = paste(setdiff(exports, aliases), collapse = ", "))
  expect_true("MVCens" %in% aliases)
})

test_that("unknown models fail closed", {
  truth <- small_truth()
  expect_error(
    mv_random("UNKNOWN", 2, truth$M, Sigma = truth$Sigma, Psi = truth$Psi),
    "must be one of"
  )
})

test_that("Monte Carlo sample sizes are positive integers", {
  truth <- small_truth()
  expect_error(
    mv_monte_carlo("MVN", sample_sizes = 2.5, replications = 1L,
                   truth = truth, workers = 1L, max_iter = 1L),
    "positive integer"
  )
  campaign <- mv_monte_carlo(
    "MVN", sample_sizes = 1L, replications = 1L,
    truth = truth, workers = 1L, max_iter = 1L
  )
  expect_equal(campaign$results$n, 1L)
})

test_that("MVREN uses unit rates as its sole identifiability convention", {
  truth <- small_truth(skew = TRUE)
  generated <- mv_random(
    "MVREN", n = 2L, M = truth$M, A = truth$A,
    Sigma = truth$Sigma, Psi = truth$Psi
  )
  expect_error(
    mv_random(
      "MVREN", n = 2L, M = truth$M, A = truth$A,
      Sigma = truth$Sigma, Psi = truth$Psi, lambda = c(2, 1)
    ),
    "every lambda_i to equal 1"
  )
  expect_error(
    mv_fit("MVREN", generated, lambda_mode = "unit_A_rows"),
    "unused argument|lambda_mode"
  )
})
