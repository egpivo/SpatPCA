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
    checkNewLocationsForSpatpcaObject(NULL, NULL),
    cat("Invalid object! Please enter a `spatpca` object")
  )
  expect_error(
    checkNewLocationsForSpatpcaObject(cv_1D, NULL),
    cat("New locations cannot be NULL")
  )
  expect_error(
    checkNewLocationsForSpatpcaObject(cv_1D, matrix(c(1, 2), ncol = 2)),
    cat("Inconsistent dimension of locations - original dimension is 1")
  )
  expect_null(checkNewLocationsForSpatpcaObject(cv_1D, x_1Dnew))
})

# Test invalid input
test_that("check input of spatpca", {
  expect_error(
    checkInputData(x = as.matrix(1), Y = Y_1D),
    cat("The number of rows of x should be equal to the number of columns of Y.")
  )
  expect_error(
    checkInputData(x = matrix(1:10, ncol = 10), Y = matrix(1)),
    cat("Number of locations must be larger than 2.")
  )
  expect_error(
    checkInputData(x = matrix(1:40, ncol = 10), Y = Y_1D),
    cat("Dimension of locations must be less than 4.")
  )
  expect_error(
    checkInputData(x = x_1D, Y = Y_1D, M = 1000),
    cat("Number of folds must be less than sample size.")
  )
  expect_null(setNumberEigenfunctions(NULL))
  expect_warning(setNumberEigenfunctions(300, 5, 100, 10))
  expect_warning(setNumberEigenfunctions(3, 5, 100, 10), NA)
})

# Test detrend
test_that("check detrending", {
  expect_equal(detrend(Y_1D, FALSE), Y_1D)
  expect_lte(sum(detrend(Y_1D, TRUE)), tol)
})
