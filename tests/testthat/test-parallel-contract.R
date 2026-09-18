test_that("independent task runner preserves task order", {
  tasks <- as.list(1:4)
  serial <- MVCens:::run_independent_tasks(tasks, identity, seed = 1, workers = 1)
  parallel <- MVCens:::run_independent_tasks(tasks, identity, seed = 1, workers = 2)
  expect_identical(serial, parallel)
})

test_that("random streams do not depend on worker count", {
  tasks <- as.list(1:6)
  worker <- function(task) stats::rnorm(3L, mean = task)
  serial <- MVCens:::run_independent_tasks(tasks, worker, seed = 871, workers = 1)
  parallel <- MVCens:::run_independent_tasks(tasks, worker, seed = 871, workers = 2)
  expect_identical(serial, parallel)
})

test_that("MVN Monte Carlo is reproducible across worker counts", {
  truth <- small_truth()
  serial <- mv_monte_carlo(
    "MVN", sample_sizes = 5L, replications = 2L, truth = truth,
    workers = 1L, max_iter = 1L, seed = 99L
  )
  parallel <- mv_monte_carlo(
    "MVN", sample_sizes = 5L, replications = 2L, truth = truth,
    workers = 2L, max_iter = 1L, seed = 99L
  )
  expect_identical(serial, parallel)
})
