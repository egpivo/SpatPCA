test_that(".onAttach sets RCPP_PARALLEL_BACKEND and displays welcome message", {
  # Save the current value of RCPP_PARALLEL_BACKEND
  original_backend <- Sys.getenv("RCPP_PARALLEL_BACKEND")
  
  # Capture package startup message
  output <- capture_message(.onAttach("SpatPCA", "SpatPCA"))
  
  # Check if RCPP_PARALLEL_BACKEND is set to "tinythread"
  expect_equal(Sys.getenv("RCPP_PARALLEL_BACKEND"), "tinythread")
  
  # Check if the welcome message is displayed
  expect_match(output$message, "Welcome to SpatPCA")
  
  # Restore the original value of RCPP_PARALLEL_BACKEND
  Sys.setenv(RCPP_PARALLEL_BACKEND = original_backend)
})

test_that(".onAttach works without interfering with CRAN submission", {
  # Simulate CRAN submission by setting NOT_CRAN to true
  Sys.setenv(NOT_CRAN = "true")
  
  # Save the current value of RCPP_PARALLEL_BACKEND
  original_backend <- Sys.getenv("RCPP_PARALLEL_BACKEND")
  
  # Call the .onAttach function
  .onAttach("SpatPCA", "SpatPCA")
  
  # Check if RCPP_PARALLEL_BACKEND is not set when NOT_CRAN is true
  expect_equal(Sys.getenv("RCPP_PARALLEL_BACKEND"), original_backend)
  
  # Reset NOT_CRAN to avoid interference with other tests
  Sys.unsetenv("NOT_CRAN")
})

test_that(".onAttach sets RCPP_PARALLEL_BACKEND for parallel processing", {
  # Save the current value of RCPP_PARALLEL_BACKEND
  original_backend <- Sys.getenv("RCPP_PARALLEL_BACKEND")
  
  # Simulate parallel processing by setting NOT_CRAN to false
  Sys.setenv(NOT_CRAN = "false")
  
  # Call the .onAttach function
  .onAttach("SpatPCA", "SpatPCA")
  
  # Check if RCPP_PARALLEL_BACKEND is set to "tinythread" for parallel processing
  expect_equal(Sys.getenv("RCPP_PARALLEL_BACKEND"), "tinythread")
  
  # Restore the original value of RCPP_PARALLEL_BACKEND
  Sys.setenv(RCPP_PARALLEL_BACKEND = original_backend)
  
  # Reset NOT_CRAN to avoid interference with other tests
  Sys.unsetenv("NOT_CRAN")
})
