.onAttach <- function(libname, pkgname) {
  if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
  }
  packageStartupMessage("Welcome to SpatPCA")
}
