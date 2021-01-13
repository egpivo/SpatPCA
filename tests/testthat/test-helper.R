defaultNumber <- RcppParallel::defaultNumThreads()
test_that("The number of cores for RcppParallel", {
  expect_error(set_cores("test"),
               "Please enter valid type - but got character")
  expect_error(set_cores(0),
               "The number of cores is not greater than 1 - but got 0")
  expect_true(set_cores(3))
  expect_true(set_cores(defaultNumber))
  expect_null(set_cores())
})
