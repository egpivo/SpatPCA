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

cv_1D <- spatpca(x = x_1D, Y = Y_1D, num_cores = num_cores)
cv_1D_fixed_K_multiple_tau2 <- spatpca(
    x = x_1D,
    Y = Y_1D,
    K = 1,
    tau2 = c(0, 1),
    num_cores = num_cores
  )
cv_1D_fixed_K_multiple_gamma <- spatpca(
    x = x_1D,
    Y = Y_1D,
    K = 1,
    gamma = c(0, 1),
    num_cores = num_cores
  )
cv_1D_fixed_K_fixed_tau1_fixed_tau2 <- spatpca(
    x = x_1D,
    Y = Y_1D,
    K = 1,
    tau1 = 10,
    tau2 = 100,
    num_cores = num_cores
  )
cv_1D_fixed_K_fixed_tau1_fixed_tau2_multiple_gamma <- spatpca(
    x = x_1D,
    Y = Y_1D,
    K = 1,
    tau1 = 0,
    tau2 = 0,
    gamma = c(0, 0.5, 1),
    num_cores = num_cores
  )
cv_1D_fixed_K_fixed_tau1 <- spatpca(
    x = x_1D,
    Y = Y_1D,
    K = 1,
    tau1 = 10,
    num_cores = num_cores
  )
estimated_eigenvalue_large_gamma <-
  spatialPrediction(
    cv_1D_fixed_K_fixed_tau1$eigenfn,
    cv_1D_fixed_K_fixed_tau1$detrended_Y,
    1000,
    cv_1D_fixed_K_fixed_tau1$eigenfn
  )$eigenvalue

estimated_signle_eigenvalue_medium_gamma <-
  spatialPrediction(
    cv_1D_fixed_K_fixed_tau1$eigenfn,
    cv_1D_fixed_K_fixed_tau1$detrended_Y,
    5,
    cv_1D_fixed_K_fixed_tau1$eigenfn
  )$eigenvalue

cv_1D_fixed_K_zero_tau1_zero_tau2 <- spatpca(
    x = x_1D,
    Y = Y_1D,
    K = 5,
    tau1 = 0,
    tau2 = 0,
    num_cores = num_cores
  )
estimated_eigenvalue_small_gamma <- spatialPrediction(
    cv_1D_fixed_K_zero_tau1_zero_tau2$eigenfn,
    cv_1D_fixed_K_zero_tau1_zero_tau2$detrended_Y,
    0.5,
    cv_1D_fixed_K_zero_tau1_zero_tau2$eigenfn
  )$eigenvalue

estimated_eigenvalue_medium_gamma <- spatialPrediction(
  cv_1D_fixed_K_zero_tau1_zero_tau2$eigenfn,
  cv_1D_fixed_K_zero_tau1_zero_tau2$detrended_Y,
  10,
  cv_1D_fixed_K_zero_tau1_zero_tau2$eigenfn
)$eigenvalue


used_number_cores <-
  as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", ""))
expected_selected_tau1_R_3.6_higher <- 0.00046416
expected_selected_tau1_R_3.6_lower <- 0.01
expected_selected_gamma_R_3.6_higher <- 0.44503397
expected_selected_gamma_R_3.6_lower <- 0.4737518

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
      abs(
        cv_1D$selected_gamma - expected_selected_gamma_R_3.6_higher
      ),
      abs(cv_1D$selected_gamma - expected_selected_gamma_R_3.6_lower)
    ),
    tol
  )
  expect_null(cv_1D$selected_K)
  expect_equal(cv_1D_fixed_K_multiple_tau2$selected_K, 1)
  expect_lte(
    abs(cv_1D_fixed_K_multiple_tau2$selected_tau1 - 0.002154435),
    tol
  )
  expect_equal(cv_1D_fixed_K_multiple_tau2$selected_tau2, 1)
  expect_equal(cv_1D_fixed_K_multiple_gamma$selected_gamma, 0)
  expect_equal(cv_1D_fixed_K_fixed_tau1_fixed_tau2$selected_tau1, 10)
  expect_equal(cv_1D_fixed_K_fixed_tau1_fixed_tau2$selected_tau2, 100)
  expect_equal(
    cv_1D_fixed_K_fixed_tau1_fixed_tau2_multiple_gamma$selected_gamma,
    0
  )
  expect_equal(cv_1D_fixed_K_fixed_tau1$selected_tau1, 10)
  expect_equal(sum(cv_1D_fixed_K_fixed_tau1$cv_score_tau1), 0)
  expect_equal(sum(estimated_eigenvalue_large_gamma), 0)
  expect_equal(sum(estimated_eigenvalue_medium_gamma), 0)
  expect_lte(abs(sum(estimated_eigenvalue_small_gamma) - 9.397599), tol)
  expect_lte(sum(estimated_signle_eigenvalue_medium_gamma) - 0.1066821, tol)
})


# Test environment settings
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
dominant_pattern_on_new_sites <-
  predictEigenfunction(cv_1D, x_new = x_1Dnew)

test_that("prediction", {
  expect_equal(ncol(prediction), 4)
  expect_equal(nrow(dominant_pattern_on_new_sites), 4)
})


# Test auxiliary function - CV with selecting K
set.seed(1234)
M <- 3
shuffle_split <- sample(rep(1:M, length.out = nrow(Y_1D)))
tau1 <- setTau1(NULL, M)
tau2 <- setTau2(NULL, M)
l2 <- setL2(tau2)
setCores(2)
cv_with_k_seleted <-
  spatpcaCVWithSelectingK(x_1D, Y_1D, M, tau1, tau2, 1, shuffle_split, 10, 1e-04, l2)

test_that("auxiliary function for selecting K", {
  expect_equal(cv_with_k_seleted$selected_K, 1)
  expect_equal(cv_with_k_seleted$cv_result$selected_gamma, 1)
  expect_equal(cv_with_k_seleted$cv_result$selected_tau1, 1)
  expect_equal(cv_with_k_seleted$cv_result$selected_tau2, 0)
})


# 3-D
set.seed(1234)
p <- 4
x <- y <- z <- as.matrix(seq(-5, 5, length = p))
d <- expand.grid(x, y, z)
Phi_3D <-
  rowSums(exp(-d^2)) / norm(as.matrix(rowSums(exp(-d^2))), "F")
Y_3D <-
  rnorm(n = 100, sd = 3) %*% t(Phi_3D) + matrix(rnorm(n = 100 * p^3), 100, p^
    3)
cv_3D <- spatpca(x = d, Y = Y_3D, num_cores = 2)
predict <-
  eigenFunction(matrix(c(0, 0, 0), 1, 3), as.matrix(d), cv_3D$eigenfn)

test_that("3D case", {
  expect_equal(dim(cv_3D$eigenfn), c(64, 2))
  expect_lte(sum(abs(predict - c(
    0.232199, 0.007501031
  ))), tol)
})
