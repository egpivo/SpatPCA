library(testthat)

if (!testthat:::on_cran()) {
  library(SpatPCA)
  Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
  test_check("SpatPCA")
}
