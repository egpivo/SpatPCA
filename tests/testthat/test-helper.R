default_number <- RcppParallel::defaultNumThreads()
test_that("The number of cores for RcppParallel", {
  expect_error(
    setCores("test"),
    "Please enter valid type - but got character"
  )
  expect_error(
    setCores(0),
    "The number of cores is not greater than 1 - but got 0"
  )
  expect_error(
    setCores(default_number + 1),
    cat("The input number of cores is invalid - default is ", default_number)
  )
  expect_true(setCores(default_number))
  expect_null(setCores())
})

x_1D <- as.matrix(seq(-5, 5, length = 10))
x_2D <- matrix(c(1, 2), ncol = 2)
test_that("Scale locations", {
  expect_equal(sum(scaleLocation(x_1D)), 5)
  expect_equal(min(scaleLocation(x_1D)), 0)
  expect_equal(max(scaleLocation(x_1D)), 1)
  expect_equal(scaleLocation(x_2D), x_2D)
})

set.seed(1234)
tol <- 1e-6
num_cores <- 2

x_1D <- as.matrix(seq(-5, 5, length = 4))
Phi_1D <- exp(-x_1D^2) / norm(exp(-x_1D^2), "F")
Y_1D <- {
  rnorm(n = 100, sd = 3) %*% t(Phi_1D) +
    matrix(rnorm(n = 100 * 4), 100, 4)
}
cv_1D <- spatpca(x = x_1D, Y = Y_1D, num_cores = num_cores)
x_1Dnew <- as.matrix(seq(6, 7, length = 4))

# Test `predict`
test_that("check new locations for a spatpca object", {
  expect_error(
    check_new_locations_for_spatpca_object(cv_1D, NULL),
    cat("New locations cannot be NULL")
  )
  expect_error(
    check_new_locations_for_spatpca_object(cv_1D, matrix(c(1, 2), ncol = 2)),
    cat("Inconsistent dimension of locations - original dimension is 1")
  )
  expect_null(check_new_locations_for_spatpca_object(cv_1D, x_1Dnew))
})

