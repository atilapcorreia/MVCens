test_that("positive-definite projection is symmetric and positive", {
  X <- matrix(c(1, 2, 2, 1), 2, 2)
  projected <- MVCens:::matrix_make_posdef(X, epsilon = 1e-6)
  expect_equal(projected, t(projected), tolerance = 1e-12)
  expect_positive_definite(projected, tolerance = 0)
})

test_that("Psi normalization enforces determinant one", {
  Psi <- matrix(c(2, 0.3, 0.3, 1.4), 2, 2)
  normalized <- MVCens:::matrix_normalize_Psi(Psi)
  expect_equal(det(normalized), 1, tolerance = 1e-10)
})

test_that("BIC helper implements the canonical formula", {
  expect_equal(MVCens:::model_bic(-10, 4, 100), 20 + 4 * log(100))
})

test_that("shared vectorization preserves the former column order", {
  X <- matrix(1:6, 2, 3)
  expect_identical(MVCens:::matrix_vectorize(X), as.vector(matrixNormal::vec(X)))
})

test_that("shared positive-definite checks reject invalid matrices", {
  expect_true(MVCens:::matrix_is_posdef(diag(2)))
  expect_false(MVCens:::matrix_is_posdef(matrix(c(1, 2, 2, 1), 2, 2)))
  expect_error(MVCens:::matrix_make_posdef(matrix(1:6, 2, 3)), "square")
})
