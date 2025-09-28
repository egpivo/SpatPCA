test_that(".onAttach only emits welcome message", {
  original_backend <- Sys.getenv("RCPP_PARALLEL_BACKEND", unset = "")
  output <- capture_message(.onAttach("SpatPCA", "SpatPCA"))
  expect_match(output$message, "Welcome to SpatPCA")
  expect_equal(Sys.getenv("RCPP_PARALLEL_BACKEND", unset = ""), original_backend)
})

test_that(".onAttach does not depend on NOT_CRAN", {
  Sys.setenv(NOT_CRAN = "true")
  original_backend <- Sys.getenv("RCPP_PARALLEL_BACKEND", unset = "")
  output <- capture_message(.onAttach("SpatPCA", "SpatPCA"))
  expect_match(output$message, "Welcome to SpatPCA")
  expect_equal(Sys.getenv("RCPP_PARALLEL_BACKEND", unset = ""), original_backend)
  Sys.unsetenv("NOT_CRAN")
})
