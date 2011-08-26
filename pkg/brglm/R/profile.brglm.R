print.profile.brglm <- function (x, ...)
{
    cat("'level' was set to", attr(x, "level"), "\n")
    cat("Methods that apply:\n")
    cat("'confint'  'plot' 'pairs'\n")
}

profile.brglm <-   function (fm,
            gridsize = 20,
            stdn = 5,
            stepsize = 0.5,
            level = 0.95,
            which = 1:length(coef(fm)),
            verbose = TRUE,
            zero.bound = 1e-08,
            scale = FALSE, ...)
{
  notNA <- !is.na(fm$coefficients)
  if (level <= 0 | level >= 1)
    stop("invalid 'level'.")
  fmML <- fitted <- update(fm,
                    method = "glm.fit",
                    control = glm.control(),
                    dispTrans = NULL,
                    correction = NULL)
  if (fm$dispTrans != "inverse") {
    if (verbose) cat("Refitting model to estimate inverse dispersion...\n")
    fm <- update(fm, dispTrans = "inverse")
  }
  Xmat <- model.matrix(fm)[, notNA]
  if (verbose)
    cat("Profiling the ordinary deviance for the corresponding ML fit...\n")
  browser()
  res1 <- profileModel(fitted = fitted,
                       gridsize = gridsize,
                       stdn = stdn,
                       stepsize = stepsize,
                       grid.bounds = NULL,
                       quantile = qchisq(level, 1),
                       objective = "ordinaryDeviance",
                       agreement = TRUE,
                       verbose = FALSE,
                       trace.prelim = FALSE,
                       which = which,
                       profTraces = TRUE,
                       zero.bound = zero.bound,
                       scale = scale,
                       dispersion = summary(fitted)$dispersion)
  fitted <- fm
  if (verbose)
    cat("Profiling the modified score statistic for the supplied fit...\n")
  res2 <- profileModel(fitted = fitted,
                       gridsize = gridsize,
                       stdn = stdn,
                       stepsize = stepsize,
                       grid.bounds = NULL,
                       quantile = qchisq(level, 1),
                       objective = "adjustedScoreStatistic",
                       agreement = TRUE,
                       verbose = FALSE,
                       trace.prelim = FALSE,
                       which = which,
                       profTraces = TRUE,
                       zero.bound = zero.bound,
                       scale = scale,
                       X = model.matrix(fitted)[, !is.na(fitted$coefficients)],
                       dispersion = fitted$dispersion)
  res <- list(profilesML = res1, profilesBR = res2)
  attr(res, "level") <- level
  class(res) <- "profile.brglm"
  res
}



