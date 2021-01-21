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
