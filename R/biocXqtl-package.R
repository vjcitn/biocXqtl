#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib biocXqtl, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

# Set before any OpenMP-linked library (e.g. torch via tensorqtl) is imported,
# to avoid "libomp already initialized" abort on macOS.
.onLoad = function(libname, pkgname) {
  Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
  Sys.setenv(OMP_NUM_THREADS = "1")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
}
