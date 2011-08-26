confint.brglm <- function (fm,
                            parm = 1:length(coef(fm)),
                            level = 0.95,
                            verbose = TRUE,
                            endpoint.tolerance = 0.001,
                            max.zoom = 100,
                            zero.bound = 1e-08,
                            stepsize = 0.5,
                            stdn = 5,
                            gridsize = 10,
                            scale = FALSE,
                            method = "smooth",
                            ci.method = "union",
                            n.interpolations = 100,
                            ...)
{
    prof <- profile.brglm(fm,
                           gridsize = 10,
                           stdn = stdn,
                           stepsize = stepsize,
                           grid.bounds = NULL,
                           level = level,
                           which = parm,
                           verbose = verbose,
                           zero.bound = zero.bound,
                           scale = scale)
    ci <- confint.profile.brglm(prof,
                                 method = method,
                                 ci.method = ci.method,
                                 endpoint.tolerance = endpoint.tolerance,
                                 max.zoom = max.zoom,
                                 n.interpolations = n.interpolations,
                                 verbose = verbose)
    drop(ci)
}


confint.profile.brglm <- function (fm,
                                    parm,
                                    level = 0.95,
                                    method = "smooth",
                                    ci.method = "union",
                                    endpoint.tolerance = 0.001,
                                    max.zoom = 100,
                                    n.interpolations = 100,
                                    verbose = TRUE, ...) {
  alpha <- 1 - attr(fm, "level")
  if (!(ci.method %in% c("union", "mean")))
    stop("Invalid 'ci.method'.")
  if (verbose)
    cat("Calculating confidence intervals for the ML fit using deviance profiles...\n")
  ci1 <- profConfint(fm$profilesML,
                     method = method,
                     endpoint.tolerance = endpoint.tolerance,
                     max.zoom = max.zoom,
                     n.interpolations = n.interpolations,
                     verbose = FALSE)
  fit <- fm$profilesBR$fit
  if (verbose) {
#            if (fit$pl | all(fit$family$link == "logit"))
#                cat("Calculating confidence intervals for the BR fit using penalized likelihood profiles...\n")
#            else cat("Calculating confidence intervals for the BR fit using modified score statistic profiles...\n")
    cat("Calculating confidence intervals for the BR fit using modified score statistic profiles...\n")
  }
  ci2 <- profConfint(fm$profilesBR,
                     method = method,
                     endpoint.tolerance = endpoint.tolerance,
                     max.zoom = max.zoom,
                     n.interpolations = n.interpolations,
                     verbose = FALSE)
  ci <- switch(ci.method,
               union = cbind(pmin(ci1[, 1],
                 ci2[, 1]), pmax(ci1[, 2], ci2[, 2])),
               mean = (ci1 + ci2)/2)
  profNames <- names(fm$profilesML$profiles)
  dimnames(ci) <- list(profNames, paste(c(alpha/2, 1 - alpha/2) *
                                        100, "%"))
  attr(ci, "profileModel object") <- NULL
  ci
}
