library(testthat)

if (!testthat:::on_cran()) {
  library(SpatPCA)
  test_check("SpatPCA")
}
