# generate 1-D data with a given seed
set.seed(1234)
tol <- 1e-6
num_cores <- 2

x_1D <- as.matrix(seq(-5, 5, length = 10))
Phi_1D <- exp(-x_1D^2) / norm(exp(-x_1D^2), "F")
Y_1D <- {
  rnorm(n = 100, sd = 3) %*% t(Phi_1D) +
    matrix(rnorm(n = 100 * 10), 100, 10)
}
# Test invalid input
test_that("check input of spatpca", {
  expect_error(
    spatpca(x = as.matrix(1), Y = Y_1D, num_cores = num_cores),
    cat("The number of rows of x should be equal to the number of columns of Y.")
  )
  expect_error(
    spatpca(x = matrix(1:10, ncol = 10), Y = matrix(1), num_cores = num_cores),
    cat("Number of locations must be larger than 2.")
  )
  expect_error(
    spatpca(x = matrix(1:10, ncol = 10), Y = Y_1D, num_cores = num_cores),
    cat("Dimension of locations must be less than 4.")
  )
  expect_error(
    spatpca(x = x_1D, Y = Y_1D, M = 1000, num_cores = num_cores),
    cat("Number of folds must be less than sample size.")
  )
})

cv_1D <- spatpca(x = x_1D, Y = Y_1D, num_cores = num_cores)

used_number_cores <- as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", ""))
expected_selected_tau1_R_3.6_higher <- 0.0021544359
expected_selected_tau1_R_3.6_lower <- 0.0004644359
expected_selected_gamma_R_3.6_higher <- 0.2137642
expected_selected_gamma_R_3.6_lower <- 0.2762986

# Test
test_that("Selected tuning parameters", {
  expect_lte(
    min(
      abs(cv_1D$selected_tau1 - expected_selected_tau1_R_3.6_higher),
      abs(cv_1D$selected_tau1 - expected_selected_tau1_R_3.6_lower)
    ),
    tol
  )
  expect_lte(abs(cv_1D$selected_tau2 - 0), tol)
  expect_lte(
    min(
      abs(cv_1D$selected_gamma - expected_selected_gamma_R_3.6_higher),
      abs(cv_1D$selected_gamma - expected_selected_gamma_R_3.6_lower)
    ),
    tol
  )
  expect_null(cv_1D$selected_K)
})

test_that("Number of threads", {
  expect_equal(num_cores, used_number_cores)
})

test_that("cross-validation plot", {
  expect_error(
    plot.spatpca("test"),
    cat("Invalid object! Please enter a `spatpca` object")
  )
  expect_true("ggplot" %in% class(plot.spatpca(cv_1D)))
})

# Test `predict`
x_1Dnew <- as.matrix(seq(6, 7, length = 4))
prediction <- predict(cv_1D, x_new = x_1Dnew)
dominant_pattern_on_new_sites <- predictEigenfunction(cv_1D, x_new = x_1Dnew)

test_that("prediction", {
  expect_equal(ncol(prediction), 4)
  expect_equal(nrow(dominant_pattern_on_new_sites), 4)
})
