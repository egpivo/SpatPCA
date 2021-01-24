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
