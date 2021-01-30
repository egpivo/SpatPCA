tol <- 1e-6
pesudo_sequence <- seq(-5, 5, length = 2)
two_dim_location <-
  as.matrix(expand.grid(x = pesudo_sequence, y = pesudo_sequence))

three_dim_location <- 
  as.matrix(expand.grid(x = pesudo_sequence, y = pesudo_sequence, z = pesudo_sequence))

# thinPlateMatrix
thin_plate_matrix_2D <- thinPlateSplineMatrix(two_dim_location)
thin_plate_matrix_3D <- thinPlateSplineMatrix(three_dim_location)
test_that("Thin-Plate Spline Matrix", {
  expect_lte(norm(thin_plate_matrix_2D, "F") - 0.362588, tol)
  expect_lte(norm(thin_plate_matrix_3D, "F") - 8.191034, tol)
})

# spatialPrediction
Phi <- matrix(c(1, 0, 0, 0), nrow = 4, ncol = 1)
new_location <- matrix(c(0.1, 0.2), nrow = 1, ncol = 2)
test_that("Eigen-function", {
  expect_lte(
    eigenFunction(new_location, two_dim_location, Phi) - 0.2352884,
    tol
  )
})
