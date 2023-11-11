default_number <- RcppParallel::defaultNumThreads()
test_that("The number of cores for RcppParallel", {
  expect_error(setCores("test"),
               "Please enter valid type - but got character")
  expect_error(setCores(0),
               "The number of cores is not greater than 1 - but got 0")
  expect_error(
    setCores(default_number + 1),
    cat(
      "The input number of cores is invalid - default is ",
      default_number
    )
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
Phi_1D <- exp(-x_1D ^ 2) / norm(exp(-x_1D ^ 2), "F")
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
    cat(
      "The number of rows of x should be equal to the number of columns of Y."
    )
  )
  expect_error(
    checkInputData(x = matrix(1:10, ncol = 10), Y = matrix(1)),
    cat("Number of locations must be larger than 2.")
  )
  expect_error(checkInputData(x = matrix(1:40, ncol = 10), Y = Y_1D),
               cat("Dimension of locations must be less than 4."))
  expect_error(
    checkInputData(x = x_1D, Y = Y_1D, M = 1000),
    cat("Number of folds must be less than sample size.")
  )
})

# Test detrend
test_that("check detrending", {
  expect_equal(detrend(Y_1D, FALSE), Y_1D)
  expect_lte(sum(detrend(Y_1D, TRUE)), tol)
})

# Test tuning parameters
test_that("check the number of eigenfunctons", {
  expect_equal(fetchUpperBoundNumberEigenfunctions(Y_1D, 5), 4)
  expect_null(setNumberEigenfunctions(NULL, Y_1D, 5))
  expect_warning(setNumberEigenfunctions(300, Y_1D, 5))
  expect_warning(setNumberEigenfunctions(3, Y_1D, 5), NA)
})

test_that("check turning parameter - tau1", {
  expect_equal(min(setTau1(NULL, 5)), 0)
  expect_equal(max(setTau1(NULL, 5)), 1)
  expect_lte(median(setTau1(NULL, 5)), 0.0004641589)
  expect_equal(median(setTau1(NULL, 1)), 1)
  expect_equal(setTau1(c(1, 2), 5), c(1, 2))
  expect_equal(setTau1(c(1, 2), 1), 2)
})

test_that("check turning parameter - tau2", {
  expect_equal(setTau2(NULL, 5), 0)
  expect_equal(setTau2(NULL, 1), 0)
  expect_equal(setTau2(c(1, 2), 5), c(1, 2))
  expect_equal(setTau2(c(1, 2), 1), 2)
})

test_that("check inner turning parameter based on tau2 - l2", {
  expect_equal(setL2(c(1, 2)), 1)
  expect_equal(setL2(-1), 1)
  expect_equal(max(setL2(1)), 1)
  expect_equal(min(setL2(1)), 0)
  expect_lte(median(setL2(1)) - 0.005994843, tol)
})

test_that("check turning parameter - gamma", {
  expect_equal(setGamma(c(1, 2)), c(1, 2))
  expect_lte(max(setGamma(NULL, Y_1D)) - 11.14708, tol)
  expect_equal(min(setGamma(NULL, Y_1D)), 0)
  expect_lte(median(setGamma(NULL, Y_1D)) - 0.06682497, tol)
})
