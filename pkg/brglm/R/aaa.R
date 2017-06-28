.onAttach <- function(libname, pkgname) {
  packageStartupMessage("'brglm' will gradually be superseded by 'brglm2' (https://cran.r-project.org/package=brglm2), which provides comprehensive utilities for mean and median bias reduction for all GLMs and for the detection of infinite estimates in binomial-response models")
}
