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
