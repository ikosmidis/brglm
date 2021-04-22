`separation.detection` <-
function (fit, nsteps = 30)
    {
     .Deprecated(msg = "'separation.detection' will be removed from 'brglm' at version 0.8. Comprehensive methods for the detection of infinite estimates in binomial-response models are provided by the 'detectseparation' R package (https://cran.r-project.org/package=detectseparation).")
    fit.class <- class(fit)[1]
    if (fit.class != "glm")
        stop("Only objects of class 'glm' are accepted.")
    eps <- .Machine$double.eps
    betasNames <- names(betas <- coef(fit))
    noNA <- !is.na(betas)
    stdErrors <- matrix(0, nsteps, length(betas))
    for (i in 1:nsteps) {
        suppressWarnings(temp.fit <- update(fit, control = glm.control(maxit = i,
            epsilon = eps)))
        stdErrors[i, noNA] <- summary(temp.fit)$coef[betasNames[noNA], "Std. Error"]
    }
    res <- sweep(stdErrors, 2, stdErrors[1, ], "/")
    colnames(res) <- names(coef(fit))
    res
}
