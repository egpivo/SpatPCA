# generate 1-D data with a given seed
set.seed(1234)
originalPar <- par(no.readonly = TRUE)
numCores <- RcppParallel::defaultNumThreads()

x_1D <- as.matrix(seq(-5, 5, length = 10))
Phi_1D <- exp(-x_1D ^ 2) / norm(exp(-x_1D ^ 2), "F")
Y_1D <- {
  rnorm(n = 100, sd = 3) %*% t(Phi_1D) +
    matrix(rnorm(n = 100 * 10), 100, 10)
}
cv_1D <- spatpca(x = x_1D, Y = Y_1D, plot.cv = TRUE,)

usedNumberCores <- as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", ""))
newPar <- par(no.readonly = TRUE)
expected_stau1_R_3.6_higher <- 0.0021544359
expected_stau1_R_3.6_lower <- 0.0004644359
expected_sgamma_R_3.6_higher <- 0.2137642
expected_sgamma_R_3.6_lower <- 0.2762986

# Test the result
tol <- 1e-6
test_that("Selected tuning parameters", {
  expect_lte(min(
    abs(cv_1D$stau1 - expected_stau1_R_3.6_higher),
    abs(cv_1D$stau1 - expected_stau1_R_3.6_lower)
  ),
  tol)
  expect_lte(abs(cv_1D$stau2 - 0), tol)
  expect_lte(min(
    abs(cv_1D$sgamma - expected_sgamma_R_3.6_higher),
    abs(cv_1D$sgamma - expected_sgamma_R_3.6_lower)
  ),
  tol)
  expect_null(cv_1D$Khat)
})

test_that("Envirorment setting", {
  expect_equal(originalPar, newPar)
})

test_that("Number of threads", {
  expect_equal(numCores, usedNumberCores)
})
