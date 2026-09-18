test_that("all governed models have specifications", {
  specs <- MVCens:::mvcens_model_specs()
  expect_setequal(
    names(specs),
    c("MVN", "MVNC", "MVSN", "MVSNC", "MVST", "MVRSN", "MVREN", "MVNIG", "MVVG")
  )
  expect_true(all(vapply(specs, inherits, logical(1), "mvcens_model_spec")))
})

test_that("estimable models expose scientific fit contracts", {
  specs <- MVCens:::mvcens_model_specs()
  estimable <- specs[c("MVN", "MVNC", "MVSN", "MVSNC", "MVST", "MVRSN", "MVREN")]
  expect_true(all(vapply(estimable, function(x) is.function(x$validate), logical(1))))
  expect_true(all(vapply(estimable, function(x) is.function(x$fit), logical(1))))
  expect_true(all(vapply(estimable, function(x) is.function(x$loglik), logical(1))))
})
