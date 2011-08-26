### Write profiling method --- output adjusted quadratic score statistic
### Need adjusted score statstic

### Deviance issue (glm evaluates null deviance wrongly)?

### Include safeBinaryRegression check

## Fix class to brglm

## Starting values: mustart, etastart?

DD <- function(expr,name, order = 1) {
  if(order < 1) stop("'order' must be >= 1")
  if(order == 1) D(expr,name)
  else DD(D(expr, name), name, order - 1)
}


brglm.fit <-
  function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
            mustart = NULL, offset = rep(0, nobs), family = gaussian(),
            control = list(), intercept = TRUE)
{
  control <- do.call("brglm.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y))
    rownames(y)
  else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv))
    stop("'family' argument seems not to be a valid family object",
         call. = FALSE)
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  dmu.deta <- family$dmu.deta
  d1afun <- family$d1afun
  d2afun <- family$d2afun
  d3afun <- family$d3afun
  d1TransDisp <- DD(control$Trans, "disp", order = 1)
  d2TransDisp <- DD(control$Trans, "disp", order = 2)
  unless.null <- function(x, if.null) if (is.null(x))
    if.null
    else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu <- unless.null(family$validmu, function(mu) TRUE)
  if (is.null(mustart)) {
    eval(family$initialize)
  }
  else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta))
      stop("invalid linear predictor values in empty model",
           call. = FALSE)
    mu <- linkinv(eta)
    if (!validmu(mu))
      stop("invalid fitted means in empty model", call. = FALSE)
    dev <- sum(dev.resids(y, mu, weights))
    w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  }
  else {
    coefold <- NULL
    eta <- if (!is.null(etastart))
      etastart
    else
      if (!is.null(start))
        if (length(start) != nvars)
          stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                        nvars, paste(deparse(xnames), collapse = ", ")),
               domain = NA)
        else {
          coefold <- start
          offset + as.vector(if (NCOL(x) == 1L)
                             x * start
          else x %*% start)
        }
      else
        family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
      stop("cannot find valid starting values: please specify some",
           call. = FALSE)
    boundary <- conv <- FALSE


    if (family$family == "binomial") {
      weights.adj <- weights + (!control$correction) * nvars/nobs
      y.adj <- (weights*y + (!control$correction) * 0.5 * nvars/nobs)/weights.adj
    }
    else {
      weights.adj <- weights
      y.adj <- y +
        if (family$family == "poisson") (!control$correction) * 0.5 * nvars/nobs
        else 0
    }
    warn <- getOption("warn")
    options(warn = -1)
    temp.fit <- glm.fit(x = x, y = y.adj, weights = weights.adj, start = start,
                          etastart = etastart, mustart = mustart,
                          offset = offset, family = family,
                          control = list(epsilon = control$epsilon,
                            maxit = 10000, trace = FALSE),
                          intercept = intercept)
    options(warn = warn)
    df.r <- temp.fit$df.residual




    if (family$family %in% c("poisson", "binomial")) {
      disp <- 1
      dispML <- 1
      transdisp <- eval(control$Trans)
    }
    else {
      if (df.r > 0) {
        if (any(temp.fit$weights == 0))
          warning("observations with zero weight not used for calculating dispersion")
        disp <- sum((temp.fit$weights * temp.fit$residuals^2)[temp.fit$weights >
                                                              0])/df.r
        prec <- 1/disp
        mu <- temp.fit$fitted
        devresids <- dev.resids(y, mu, weights)
        for (i in 1L:100) {
          zetas <- -weights*prec
          expDevresids <- weights*d1afun(zetas)
          precScores <- -1/2*sum(devresids - expDevresids)
          precInfo <- 1/2*sum(weights^2*d2afun(zetas))
          disp <- disp - disp^2 * precScores/precInfo
          prec <- 1/disp
          transdisp <- eval(control$Trans)
          if (critDispML <- sum(abs(precScores)) < 1e-08) {
            dispML <- disp
            break
          }
        }
        if (!critDispML)
          warning("the ML estimate of the dispersion could not be calculated. An alternative estimate had been used as starting value.")
      }
      else {
        disp <- 1 ## A random value
        dispML <- NA
        transdisp <- eval(control$Trans)
      }
    }
    nacoefs <- is.na(temp.fit$coefficients)
    coefnames <- names(temp.fit$coefficients)
    X <- x[, coefnames[!nacoefs], drop = FALSE]
    coefs <- temp.fit$coefficients[!nacoefs]
    nvarsRed <- length(coefs)
    adjScores <- rep(NA, nvars + 1)
    names(adjScores) <- c("Transformed Dispersion", coefnames)
    if (control$correction) {
      control$maxit <- 1
      control$slowit <- 1
    }
    objCur <- .Machine$integer.max
    for (iter in 1L:control$maxit) {
      halfstep <- 0
      testhalf <- TRUE
      coefsPrev <- coefs
      objPrev <- objCur
      while (testhalf & halfstep < 10) {
        good <- weights > 0
        eta <- drop(X %*% coefs + offset)
        mu <- linkinv(eta)
        dmu <- mu.eta(eta)
        ddmu <- dmu.deta(eta)
        varmu <- variance(mu)[good]
        if (any(is.na(varmu)))
          stop("NAs in V(mu)")
        if (any(varmu == 0))
          stop("0s in V(mu)")
        if (any(is.na(dmu[good])))
          stop("NAs in d(mu)/d(eta)")
        good <- (weights > 0) & (dmu != 0)
        if (all(!good)) {
          conv <- FALSE
          warning("no observations informative at iteration ",
                  iter)
          break
        }
        ## Estimation of coefficients
        w <- rep.int(0, nobs)
        w[good] <- weights[good]*dmu[good]^2/varmu
        Wx <- sqrt(w[good]) * X[good, ]
        XWXinv <- chol2inv(chol(crossprod(Wx)))
        hats <- diag(X%*%XWXinv%*%t(w * X))
        nvalid <- sum(good)
        adjExpComp <- 0.5 * hats[good] * ddmu[good]/dmu[good] * X[good, ]
        adjExp <- .Internal(colSums(adjExpComp, nvalid, nvarsRed, FALSE))
        scoresComp <- weights[good] * dmu[good]/varmu * (y - mu)[good] * X[good,]
        scores <- .Internal(colSums(scoresComp, nvalid, nvarsRed, FALSE))
        firstOrderBiasBeta <- - disp * XWXinv %*% adjExp
        adjScoresBeta <- scores/disp + adjExp
        objCur <- sqrt(disp * adjScoresBeta%*%XWXinv%*%adjScoresBeta)
#        XWX <- crossprod(Wx)
#        XWX <- XWX + min(eigen(XWX)$values)*diag(nrow(XWX))
#        coefs <- coefsPrev + 2^(-halfstep) * (drop(solve(XWX)%*%scores) +
#                                             disp * solve(XWX) %*% adjExp)
        coefs <- coefsPrev + 2^(-halfstep) * (drop(XWXinv%*%scores) -
                                              firstOrderBiasBeta)
        halfstep <- halfstep + 1
        testhalf <- objCur > objPrev
        #if (halfstep == 1) { cat("##########\n")}
        #cat("iteration:", iter,
        #    "halving factor:", halfstep,
        #    "previous:", objPrev,
        #    "current:", objCur,
        #    "test:", testhalf, "\n")
      }
      ## Dispersion estimation
      devresids <- dev.resids(y, mu, weights)
      if (family$family %in% c("poisson", "binomial")) {
        disp <- 1
        transdisp <- eval(control$Trans)
        adjScoresTransDisp <- NA
      }
      else {
        if (df.r > 0) {
          zetas <- -weights/disp
          expDevresids <- weights*d1afun(zetas)
          dispScores <- sum(devresids - expDevresids)/(2*disp^2)
          dispInfo <- (denom1 <- sum(weights^2*d2afun(zetas)))/(2*disp^4)
          adjExpDisp <- (nvarsRed - 2)/(2*disp) +
            sum(weights^3 * d3afun(zetas))/(2*disp^2*denom1)
          firstOrderBiasDisp <- - adjExpDisp/dispInfo
          transDispScores <- dispScores/(firstDeriv <- eval(d1TransDisp))
          transDispInfo <- dispInfo/firstDeriv^2
          firstOrderBiasTransDisp <- firstDeriv*firstOrderBiasDisp +
            eval(d2TransDisp)/(2*dispInfo)
          adjScoresTransDisp <- transDispScores -
            transDispInfo * firstOrderBiasTransDisp
          transdisp <- transdisp + control$slowit *
            (transDispScores/transDispInfo - firstOrderBiasTransDisp)
          disp <- eval(control$inverseTrans)
        }
        else {
          disp <- 1 # No effect to the adjusted scores
          transdisp <- eval(control$Trans)
          adjScoresTransDisp <- NA
        }
      }
      adjScores[c(TRUE, !nacoefs)] <- c(adjScoresTransDisp, adjScoresBeta)
      if (control$trace) {
        cat("Iteration:", iter, "\n")
        cat("Adjusted scores:", adjScores, "\n")
        browser()
      }
      if (conv <- sum(abs(adjScores), na.rm = TRUE) < control$epsilon) {
        break
      }
    }
    if (df.r == 0) disp <- NaN
    coef <- rep(NA, length(temp.fit$coefficients))
    coef[!nacoefs] <- coefs
    ## More associated with
    if ((!conv) & (!control$correction))
      warning("brglm.fit: algorithm did not converge", call. = FALSE)
    ## Boundary?
    if (boundary)
      warning("brglm.fit: algorithm stopped at boundary value",
              call. = FALSE)
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
      if (any(mu > 1 - eps) || any(mu < eps))
        warning("brglm.fit: fitted probabilities numerically 0 or 1 occurred",
                call. = FALSE)
    }
    if (family$family == "poisson") {
      if (any(mu < eps))
        warning("brglm.fit: fitted rates numerically 0 occurred",
                call. = FALSE)
    }
    Qr <- qr(sqrt(w)*x)
    temp.fit$qr <- as.matrix(Qr$qr)
    temp.fit$pivot <- Qr$pivot
    temp.fit$rank <- Qr$rank
    temp.fit$qraux <- Qr$qraux
    xxnames <- xnames[temp.fit$pivot]
    residuals <- (y - mu)/mu.eta(eta)
    nr <- min(nvalid, nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1L:nr, 1L:nvars] <- temp.fit$qr[1L:nr, 1L:nvars]
    }
    else Rmat <- temp.fit$qr[1L:nvars, 1L:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    colnames(temp.fit$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  wt <- rep.int(0, nobs)
  wt[good] <- w
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  if (!EMPTY)
    names(temp.fit$effects) <- c(xxnames[seq_len(temp.fit$rank)], rep.int("",
                                                                sum(good) - temp.fit$rank))
  wtdmu <- if (intercept)
    sum(weights * y)/sum(weights)
  else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- if (EMPTY)
    0
  else temp.fit$rank
  resdf <- n.ok - rank
  temp.fit$deviance <- sum(devresids)
  aic.model <- aic(y, n, mu, weights, temp.fit$deviance) + 2 * rank
  list(coefficients = coef,
       residuals = residuals,
       fitted.values = mu,
       effects = if (!EMPTY) temp.fit$effects,
       R = if (!EMPTY) Rmat,
       rank = rank,
       qr = if (!EMPTY) structure(temp.fit[c("qr", "rank",
         "qraux", "pivot", "tol")], class = "qr"),
       family = family,
       linear.predictors = eta,
       deviance = temp.fit$deviance,
       aic = aic.model,
       null.deviance = nulldev,
       iter = iter,
       weights = wt,
       prior.weights = weights,
       df.residual = resdf,
       df.null = nulldf,
       y = y,
       converged = conv,
       boundary = boundary,
       dispersion = disp,
       dispersionML = dispML,
       transDispersion = transdisp,
       adjustedScores = adjScores,
       dispTrans = control$dispTrans,
       cov.unscaled = XWXinv,
       class = "brglm")
}

summary.brglm <- function (object, dispersion = object$dispersion,
                            correlation = FALSE, symbolic.cor = FALSE,
                            ...) {
  summary.glm(object, dispersion = dispersion,
              correlation = correlation,
              symbolic.cor = symbolic.cor, ...)
}

customTrans <- list(Trans = expression(disp),
                    inverseTrans = expression(transdisp))

