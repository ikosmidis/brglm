.onAttach <- function(libname, pkgname) {
  packageStartupMessage("'brglm' will gradually be superseded by the 'brglm2' R package (https://cran.r-project.org/package=brglm2), which provides utilities for mean and median bias reduction for all GLMs.\n Methods for the detection of separation and infinite estimates in binomial-response models are provided by the 'detectseparation' R package (https://cran.r-project.org/package=detectseparation).")
}
