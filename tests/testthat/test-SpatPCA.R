# generate 1-D data with a given seed
set.seed(1234)
x_1D <- as.matrix(seq(-5, 5, length = 10))
Phi_1D <- exp(-x_1D ^ 2) / norm(exp(-x_1D ^ 2), "F")
Y_1D <- {
  rnorm(n = 100, sd = 3) %*% t(Phi_1D) +
    matrix(rnorm(n = 100 * 10), 100, 10)
}
cv_1D <- spatpca(x = x_1D, Y = Y_1D)
expected_stau1_R_3.6_higher <- 0.0021544359
expected_stau1_R_3.6_lower <- 0.0004644359

expected_sgamma_R_3.6_higher <- 0.2137642
expected_sgamma_R_3.6_lower <- 0.2762299

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
