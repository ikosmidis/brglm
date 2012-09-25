## Starting values is essential

### Clean up calls to names in bpolr

## Fix the odd handling of nominal effects, ie. sign...

## .dat is redundant object in bpolr

## Check thisCall$offset

## Location, Scale, Nominal
## organise <- function(formula, data, weights, subset, na.action,
##                      contrasts = NULL) {
organise <- function(formula, contrasts = NULL, mf,
#                     offsetLinear = NULL,
#                     offsetScale = NULL,
                     drop.only.mid.categories = FALSE,
                     drop.empty.categories = FALSE) {
  require(gnm)
  require(Formula)

  ## To be used for dat resulted from expandCategorical
  groupDat <- function(dat, response, weights) {
    response1 <- dat[, response]
    datNo <- dat[,-match(c(weights), names(dat))]
    datNo[[response]] <- NULL
    datNo[[response]] <- response1
    ord <- do.call("order", datNo)
    datord <- datNo[ord,]
    dupp <- duplicated(datord)
    inds <- which(!dupp)
    reps <- c(diff(inds), nrow(datord) - rev(inds)[1] + 1)
    inds1 <- rep(1:length(reps), reps)
    wt <- c(tapply(dat[[weights]][ord], inds1, sum))
    datt <- datord[inds,]
    datt[[weights]] <- wt
    datt
  }

  cleanAliasing <- function(mat) {
    matQR <- qr(mat, tol = 1e-08, LAPACK = FALSE)
    R <- matQR$rank
    aliased <- logical(ncol(mat))
    names <- colnames(mat)
    if (R != ncol(mat)) {
      mat <- mat[, matQR$pivot[1:matQR$rank], drop = FALSE]
      aliased[-matQR$pivot[1:matQR$rank]] <- TRUE
      if (R != qr(mat)$rank) {
        stop("Failed to find a full rank design matrix", call. = FALSE)
      }
    }
    attr(mat, "aliased") <- structure(aliased, .Names = names)
    attr(mat, "variables") <- names
    mat
  }

  thisCall <- match.call()
  ## if (missing(formula)) stop("A specification of location is absent")
  ## if (missing(data)) {
  ##   data <- environment(formula)
  ## }
  ## mf <- match.call(expand.dots = FALSE)
  ## mf$contrasts <- NULL
  ## mf$drop.unused.levels <- TRUE
  ## oformula <- as.formula(formula)
  ## formula <- as.Formula(formula)
  ## environment(formula) <- parent.frame()
  ## if (length(formula)[2L] < 2L) {
  ##   formula <- as.Formula(formula(formula), ~ 1, ~ 1)
  ##   simple_formula <- TRUE
  ## }
  ## if (length(formula)[2L] == 2) {
  ##   formula <- as.Formula(formula(formula), ~ 1)
  ##   simple_formula <- TRUE
  ## }

  ## else {
  ##   if (length(formula)[2L] > 3L) {
  ##     formula <- Formula(formula(formula, rhs = 1:3))
  ##     warning("formula must not have more than three RHS parts")
  ##   }
  ##   simple_formula <- FALSE
  ## }
  ## mf$formula <- formula
  ## mf[[1L]] <- as.name("model.frame")
  ## mf <- eval(mf, parent.frame())

  .dataCurrent <- mf
  termsmf <- attr(.dataCurrent, "terms")
  nam <- names(.dataCurrent)
  responseName <- names(model.part(formula, data = .dataCurrent, lhs = 1L))

  if (!("(weights)" %in% nam)) {
    .dataCurrent[["(weights)"]] <- rep(1, nrow(.dataCurrent))
    nam <- c(nam, "(weights)")
  }

  if (!is.factor(model.response(.dataCurrent))) {
    stop("The response needs to be a factor")
  }

  mtL <- terms(formula, data = .dataCurrent, rhs = 1L)
  interceptOnly <- identical(all.vars(mtL), responseName)
  if (interceptOnly) {
    .dataCurrent$TEMPORARYVARIABLE <- 1
    nam <- c(nam, "TEMPORARYVARIABLE")
  }


  ## ### Offsets
  ## .dataCurrent$offsetLinear <- offsetLinear
  ## .dataCurrent$offsetScale <- offsetScale


  .dataCurrent <- expandCategorical(.dataCurrent, responseName, group = FALSE)
  .dataCurrent[["(weights)"]] <- c(.dataCurrent[["(weights)"]] * .dataCurrent[["count"]])
  inds <- match(c(nam[-1], nam[1]), names(.dataCurrent))
  .dataCurrent <- groupDat(.dataCurrent[, inds], responseName, "(weights)")
  .dataCurrent <- .dataCurrent[nam]

  if (interceptOnly) {
    .dataCurrent$TEMPORARYVARIABLE <- NULL
  }

  ## Remove any non-observed categories
  attr(.dataCurrent, "terms") <- termsmf

  if (drop.empty.categories) {
    empty <- tapply(model.weights(.dataCurrent),
                    Y <- model.response(.dataCurrent),
                    sum) == 0
    ## If you need to drop all but one category then kepp the
    ## non-empty one with one and only one empty category
    nlev <- nlevels(Y)
    if (any(empty)) {
      if (drop.only.mid.categories & (empty[1] | empty[nlev])) {
        empty[1] <- empty[nlev] <- FALSE
      }
      if (sum(!empty) < 2) {
        empty[max(which(empty))] <- FALSE
      }
      empty <- which(empty)
      .dataCurrent <- .dataCurrent[!match(Y,  levels(Y)[empty], nomatch = 0), ]
      .dataCurrent <- droplevels(.dataCurrent,
                                 except = !match(names(.dataCurrent), responseName, nomatch = 0))
      attr(.dataCurrent, "terms") <- termsmf
    }
  }

  ## Get model terms, responses and so on
  mt <- terms(formula, data = .dataCurrent)
  mtL <- terms(formula, data = .dataCurrent, rhs = 1L)
  mtS <- delete.response(terms(formula, data = .dataCurrent, rhs = 2L))
  mtN <- delete.response(terms(formula, data = .dataCurrent, rhs = 3L))
  Y <- model.response(.dataCurrent, "any")
  N <- length(Y)
  if (N < 1)
    stop("empty model")

  ## Useful integers
  lev <- levels(Y)
  nlev <- length(lev)
  q <- nlev - 1
  NcovClass <- N/nlev
  Nobs <- NcovClass*q
  rownames(.dataCurrent) <- paste(rep(1:NcovClass, each = nlev),
                                  rep(1:nlev, times = NcovClass), sep = ".")

  ## Set-up model matrices
  XL <- model.matrix(mtL, .dataCurrent, contrasts = contrasts)
  XN <- model.matrix(mtN, .dataCurrent, contrasts = contrasts)
  XS <- model.matrix(mtS, .dataCurrent, contrasts = contrasts)

  XLInt <-  match("(Intercept)", colnames(XL), nomatch = 0)
  XNInt <-  match("(Intercept)", colnames(XN), nomatch = 0)
  XSInt <-  match("(Intercept)", colnames(XS), nomatch = 0)

  ## If Intercept is not in the model matrix then add it for checking
  ## the aliasing later, and also give warnings if this happens in the
  ## location specification
  if (XLInt <= 0) {
    XL <- cbind(`(Intercept)` = rep(1, nrow(XL)), XL)
    warning("an intercept is needed and assumed in the location specification")
  }
  if (XNInt <= 0) {
    XN <- cbind(`(Intercept)` = rep(1, nrow(XN)), XN)
    ## warning("an intercept is needed and assumed in the nominal specification")
  }
  if (XSInt <= 0) {
    XS <- cbind(`(Intercept)` = rep(1, nrow(XS)), XS)
    ## warning("an intercept is needed and assumed in the scale specification")
  }

  ## Take care of aliasing within nominal
  XN <- cleanAliasing(XN)
  aliasedN <- attr(XN, "aliased")

  ## Take care of nominal-location aliasing
  XNL <- cleanAliasing(cbind(XN, XL))
  aliasedNL <- attr(XNL, "aliased")
  aliasedL <- aliasedNL[-seq.int(length.out = ncol(XN))]
  XL <- XL[, !aliasedL, drop = FALSE]

  ## Take care of nominal-scale aliasing
  XNS <- cleanAliasing(cbind(XN, XS))
  aliasedNS <- attr(XNS, "aliased")
  aliasedS <- aliasedNS[-seq.int(length.out = ncol(XN))]
  colnames(XS) <- paste("scale", names(aliasedS), sep = ".")
  XS <- XS[, !aliasedS, drop = FALSE]

  ## Remove intercept from nominal (forced to be in column 1)
  XN <- XN[, -1, drop = FALSE]

  ## Set parameter names and overall aliasing vector
  .polrNam <- paste(lev[-nlev], lev[-1L], sep = "|")
  aliasedL <- aliasedL[-1]
  aliasedS <- aliasedS[-1]
  aliasedN <- aliasedN[-1]
  namesN <- names(aliasedN)
  namesNN <- NULL
  aliasedNN <- NULL
  for (i in namesN) {
    namesNN <- c(namesNN, paste(.polrNam, i, sep = "."))
    aliasedNN <- c(aliasedNN, rep(aliasedN[i], q))
  }
  aliased <- c(aliasedL, aliasedNN, logical(q), aliasedS)

  alphaNames <- .polrNam
  betaNames <- c(if (length(aliasedL)) names(aliasedL) else NULL,
                 namesNN)
  tauNames <- if (length(aliasedS)) paste("scale", names(aliasedS), sep = ".") else NULL
  names(aliased) <- c(betaNames, alphaNames, tauNames)

  ## Set-up appropriate model matrices --- include .polr
  XL <- XL[-seq(nlev, N, nlev), , drop = FALSE]
  .polr <- col(matrix(0, nrow(XL), q)) == rep(1:q, N/nlev) - 1L + 1L
  colnames(.polr) <- .polrNam

  XN <- XN[-seq(nlev, N, nlev), , drop = FALSE]
  pN <- ncol(XN)
  nominalNames <- colnames(XN)
  XNNew <- NULL
  if (pN > 0) {
      for (i in 1:pN) {
        Newmat <- XN[, i]*.polr
        colnames(Newmat) <- paste(.polrNam, nominalNames[i], sep = ".")
        XNNew <- cbind(XNNew, Newmat)
      }
  }

  XL <- cbind(XL, if (pN > 0) XNNew else NULL)
  XLinear <- cbind(-XL, .polr)

  XS <- XS[-seq(nlev, N, nlev), , drop = FALSE]

  ## Handle weights
  weights <- model.weights(.dataCurrent)
  weights <- as.vector(weights)
  names(weights) <- rownames(.dataCurrent)

  ## Handle offsets
  expand_offset <- function(offset) {
    if (is.null(offset))
      offset <- 0
    if (length(offset) == 1)
      offset <- rep.int(offset, N)
    as.vector(offset)
  }

  offsetL <- expand_offset(model.offset(model.part(formula,
                                                   data = .dataCurrent,
                                                   rhs = 1L,
                                                   terms = TRUE)))
  offsetS <- expand_offset(model.offset(model.part(formula,
                                                   data = .dataCurrent,
                                                   rhs = 2L,
                                                   terms = TRUE)))
  offsetN <- expand_offset(model.offset(model.part(formula,
                                                   data = .dataCurrent,
                                                   rhs = 3L,
                                                   terms = TRUE)))
  if (!is.null(thisCall$offset)) {
    offsetXL <- offsetXL + expand_offset(.dataCurrent[, "(offset)"])
  }

  ## offsetL will be substracted
  ## offsetS will be added
  offsetL <- offsetL[-seq(nlev, N, nlev)]
  offsetS <- offsetS[-seq(nlev, N, nlev)]
  offsetN <- offsetN[-seq(nlev, N, nlev)]

  freqs <- matrix(.dataCurrent[["(weights)"]], nrow = nlev)
  list(XL = XL,
       XLinear = XLinear,
       XS = XS,
       weights = weights,
       freqs = freqs,
       totals = colSums(freqs),
       dataReshaped = .dataCurrent,
       dataOriginal = data,
       offsetL = offsetL,
       offsetS = offsetS,
       offsetN = offsetN,
       N = N,
       q = q,
       nlev = nlev,
       NcovClass = NcovClass,
       pN = ncol(XN),
       pXL = ncol(XL),
       pXS = ncol(XS),
       aliased = aliased,
       locationFormula = formula(formula, lhs = 1L, rhs = 1L),
       scaleFormula = formula(formula, lhs = 0L, rhs = 2L),
       nominalFormula = formula(formula, lhs = 0L, rhs = 3L),
       alphaNames = alphaNames,
       betaNames = betaNames,
       tauNames = tauNames)
}



bpolr <- function(formula,
                  data,
                  weights,
                  start = NULL,
                  subset,
                  na.action,
                  contrasts = NULL,
                  model = TRUE,
                  link = c("logit", "probit", "cloglog", "cauchit"),
                  method = c("ML", "BR", "BC"),
                  maxit = 500,
                  epsilon = 1e-08,
                  history = TRUE,
                  trace = FALSE,
                  drop.empty.categories = TRUE,
                  stepFactor.init = 0,
#                  offsetLinear = NULL,
#                  offsetScale = NULL,
                  ...) {
  require(ordinal)
  require(Formula)

  if (missing(formula)) stop("A specification of location is absent")

  if (missing(data)) {
    data <- environment(formula)
  }

  link <- match.arg(link)
  pfun <- switch(link,
                 logit = make.link("logit")$linkinv,
                 probit = make.link("probit")$linkinv,
                 cloglog = make.link("cloglog")$linkinv,
                 cauchit = make.link("cauchit")$linkinv)
  dfun <- switch(link,
                 logit = make.link("logit")$mu.eta,
                 probit = make.link("probit")$mu.eta,
                 cloglog = make.link("cloglog")$mu.eta,
                 cauchit = make.link("cauchit")$mu.eta)
  ddfun <- switch(link,
                  logit = function(eta) dfun(eta)*(1 - 2*pfun(eta)),
                  probit = function(eta) -eta*dfun(eta),
                  cloglog = function(eta) dfun(eta)*(1 - exp(eta)),
                  cauchit = function(eta) -2*pi*eta*dfun(eta)^2)

  method <- match.arg(method)


  thisCall <- match.call()
  if (missing(formula)) stop("A specification of location is absent")
  if (missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  mf$contrasts <- mf$start <-
    mf$model <- mf$link <- mf$method <- mf$maxit <- mf$epsilon <-
      mf$history <- mf$trace <- mf$... <- mf$drop.empty.categories <-
        mf$stepFactor.init <- NULL
  mf$drop.unused.levels <- TRUE
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  ##  environment(formula) <- parent.frame()
  if (length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1, ~ 1)
    ## simple_formula <- TRUE
  }
  if (length(formula)[2L] == 2) {
    formula <- as.Formula(formula(formula), ~ 1)
    ## simple_formula <- TRUE
  }

  else {
    if (length(formula)[2L] > 3L) {
      formula <- Formula(formula(formula, rhs = 1:3))
      warning("formula must not have more than three RHS parts")
    }
    ## simple_formula <- FALSE
  }
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  object <- organise(formula, contrasts, mf,
#                     offsetLinear = offsetLinear,
#                     offsetScale = offsetScale,
                     drop.empty.categories = drop.empty.categories,
                     drop.only.mid.categories = (method == "BR"))

  ## subset <- if (missing(subset)) NULL else subset
  na.action <- if (missing(na.action)) NULL else na.action
  ## weights <- if (missing(weights)) NULL else weights
  ## object <- organise(formula = formula,
  ##                    data = data,
  ##                    weights = weights,
  ##                    subset = subset,
  ##                    na.action = na.action,
  ##                    contrasts = contrasts)



  ## A clm call for starting values later; used only if start is not
  ## specified and plan is to get rid of it
  Mclm <- match.call(expand.dots = FALSE)
  Mclm$start <- Mclm$method <- Mclm$link <- Mclm$model <-
    Mclm$maxit <- Mclm$epsilon <-
    Mclm$history <- Mclm$trace <- Mclm$... <- NULL

  ############################
  ## Maybe work within object?
  ############################

  ## Useful integers
  nlev <- object$nlev
  q <- object$q
  N <- object$N
  .dat <- object$dataReshaped
  NcovClass <- object$NcovClass
  Nobs <- object$Nobs
  pXS <- object$pXS
  pXL <- object$pXL

  ## Frequencies and totals
  freqs <- object$freqs
  totals <- object$totals

  ## Get model matrices
  XL <- object$XL
  XLinear <- object$XLinear
  XS <- object$XS

  ## Get offsets
  offsetL <- object$offsetL
  offsetS <- object$offsetS

  ## Aliasing vector and names
  aliased <- object$aliased
  alphaNames <- object$alphaNames
  betaNames <- object$betaNames
  tauNames <- object$tauNames

  ############################

  coefNames <- c(betaNames, alphaNames, tauNames)

  ## Necessary quantities
  inds <- sapply(1:NcovClass,
                 function(i) (1 + (i - 1)*q):(i*q))
  if (q > 1) {
    mitre <- diag(q) - cbind(0, rbind(diag(q - 1), 0))
  }
  else {
    dim(inds) <- c(1, NcovClass)
    mitre <- as.matrix(1)
  }

  fitFun <- function(pars, deriv = 1L) {
    beta <- pars[seq.int(length.out = pXL)]
    alpha <- pars[seq.int(length.out = q) + pXL]
    tau <- pars[seq.int(length.out = pXS) + q + pXL]
    etas <- matrix(c(drop(XLinear %*% c(beta, alpha)) - offsetL)/
                   exp(drop(XS %*% tau + offsetS)), nrow = q)
    Z <- cbind(XLinear/exp(drop(XS %*% tau  + offsetS)), if (pXS > 0) -XS*c(etas) else NULL)
    cumprob <- pfun(etas)
    gs <- dfun(etas)
    prob <- apply(cumprob, 2, function(x) diff(c(0, x, 1)))
    ggs <- if (deriv >= 1L) ggs <- ddfun(etas) else NULL
    list(Z = Z, etas = etas, prob = prob, gs = gs, ggs = ggs,
         beta = beta, alpha = alpha, tau = tau)
  }

  infoFun <- function(pars, fit = fitFun(pars), inverse = FALSE) {
    with(fit, {
      FisherInfo <- 0
      for (i in 1:NcovClass) {
        Zr <- Z[inds[,i], , drop = FALSE]
        ## CAN YOU AVOID THE IF BELOW?
        invCovr <- 1/prob[nlev, i] + if (q == 1) 1/prob[1, i, drop = FALSE] else diag(1/prob[1:q, i])
        Dr <- mitre * gs[, i]
        Wr <- totals[i] * Dr %*% invCovr %*% t(Dr)
        FisherInfo <- FisherInfo + t(Zr) %*% Wr %*% Zr
      }
      if (inverse) try(solve(FisherInfo), silent = TRUE) else FisherInfo
    })
  }

  gradFun <- function(pars, fit = fitFun(pars)) {
    with(fit, -colSums(c(gs * apply(freqs/prob, 2, diff)) * Z))
  }

  biasFun <- function(pars,
                   fit = fitFun(pars, deriv = 2L),
                   vcov = infoFun(pars, fit = fit, inverse = TRUE)) {
    if (inherits(vcov, "try-error")) {
      list(bias = rep(NA_real_, pXL + q + pXS),
           adjustment = rep(NA_real_, pXL + q + pXS))
    }
    else {
      adj <- 0
      with(fit, {
        cc <- matrix(0, nrow(prob), ncol(prob))
        if (pXS == 0) {
          for (i in 1:NcovClass) {
            Zr <- Z[inds[,i], , drop = FALSE]
            Ar <- Zr %*% vcov %*% t(Zr)
            cr <- 0.5 * totals[i] * diag(Ar) * ggs[, i]
            cc[, i] <- c(cr[1], diff(cr), -cr[q])
          }
        }
        else {
          invFtt <- vcov[pXL + q + 1:pXS, pXL + q + 1:pXS]
          invFat <- vcov[pXL + 1:q, pXL + q + 1:pXS]
          invFbt <- vcov[1:pXL, pXL + q + 1:pXS]
          taus <- pars[q + pXL + 1:pXS]
          for (i in 1:NcovClass) {
            Zr <- Z[inds[,i], , drop = FALSE]
            Ar <- Zr %*% vcov %*% t(Zr)
            cr <- 0.5 * totals[i] * diag(Ar) * ggs[, i]
            cc[, i] <- c(cr[1], diff(cr), -cr[q])
            etar <- etas[, i]
            invCovr <- 1/prob[nlev, i] + if (q == 1) 1/prob[1, i, drop = FALSE] else diag(1/prob[1:q, i])
            Dr <- mitre * gs[, i]
            Wr <- totals[i]*Dr%*%invCovr%*%t(Dr)
            vr <- XS[inds[1, i], ]
            xr <- XL[inds[1, i], ]
            adj <- adj + 0.5 * (vr %*% invFtt %*% vr) * c(Wr %*% etar) * Zr -
              exp(-sum(taus * vr)) * c( Wr %*% invFat %*% vr -
                                   if (pXL > 0)
                                       xr %*% invFbt %*% vr * rowSums(Wr)
                                   else 0)*Zr
          }
          adj <- colSums(adj)
        }
        adjustment <- adj -colSums(c(gs * apply(cc/prob, 2, diff)) * Z)
        bias <- -vcov %*% adjustment
        list(bias = bias, adjustment = adjustment)
      })
    }
  }

  if (is.null(start)) {
    ## Starting values
    ## Add something small to zero weights to avoid the
    ## possibility of boundary estimates
    ### Hmmm... maybe get starting values using many binomial regressions?
    www <- object$weights
    www[www == 0] <- 0.5*(pXL + q + pXS)/sum(www)
    Mclm[[1L]] <- as.name("clm")
    Mclm$formula <- object$locationFormula
    Mclm$scale <- object$scaleFormula
#    Mclm$offsetLinear <- NULL
#    Mclm$offsetScale <- NULL
    Mclm$nominal <- object$nominalFormula
    Mclm$link <- as.name("link")
    Mclm$weights <-as.name("www")
    Mclm$control <- list(maxIter = 2)
    datS <- .dat
    off <- grep("offset", names(datS))
    names(datS)[off] <- sub("\\)", "", sub("offset\\(", "", names(datS)[off]))
    datS$www <- www
    Mclm$data <- as.name("datS")
    options(warn = -1)
    clmObject <- eval(Mclm)
    options(warn = 0)
    ## zeta is what is called tau here
    pars <- c(clmObject$beta, -clmObject$alpha[-c(1:q)],
              clmObject$alpha[1:q], clmObject$zeta)
  }
  else {
    if (pXS <= 0) {
      ## browser()
      if (length(start)!=(pXL + q))
        stop("'start' is not of the correct length", call. = FALSE)
    }
    else {
      if (length(start)!=(pXL + pXS + q))
        stop("'start' is not of the correct length", call. = FALSE)
    }
    pars <- start
  }

  step <- .Machine$integer.max
  historyEstimates <- NULL
  niter <- 0
  failedInv <- FALSE

  #### Fix me --- too much names here...
  names(pars) <- coefNames
  pars <- pars[!aliased]

  parsPrev <- pars
  ## if maxit is 0 then simply evaluate everything at start
  if (maxit > 0) {
    ## Increase maxit by 1; otherwise if maxit is 1 nothing happens
    ## because we move back to parsPrev after the end of the iteration
    for (niter in 1:(maxit + 1)) {
      stepPrev <- step
      ## if stepFactor is set to zero the we have the most agressive
      ## first step, but we can be conservative here if
      ## stepFactor.init is set to a value greater than 0
      stepFactor <- stepFactor.init
      testhalf <- TRUE
      if (history) {
        historyEstimates <- cbind(historyEstimates, pars)
      }
      while (testhalf & stepFactor < 16) {
        fit <- fitFun(pars, deriv = 2L)
        scores <- gradFun(pars, fit = fit)
        infoInv <- infoFun(pars, fit = fit, inverse = TRUE)
        if (failedInv <- inherits(infoInv, "try-error")) {
          warning("failed to invert the information matrix: iteration stopped prematurely")
          break
        }
        ## If invoInv cannot be inverted then the following is not evaluated
        parsPrev <- pars
        if (method == "BR") {
          biasObject <- biasFun(pars, fit = fit, vcov = infoInv)
          bias <- biasObject$bias
          adjustment <- biasObject$adjustment
        }
        else {
          bias <- 0
          adjustment <- 0
        }
        if (trace) {
          cat("===========================\n")
          cat("Iteration: ", niter, "step factor", 2^(-stepFactor),  "\n")
          cat("Inverse condition number:",
              format(1/kappa(infoInv, exact = TRUE), nsmall = 5), "\n")
          cat("Scores:\n", format(scores + adjustment, scientific = TRUE), "\n")
          cat("===========================\n")
          browser()
        }
        pars <- pars +
          2^(-stepFactor) * (step <- infoInv %*% scores - bias)
        stepFactor <- stepFactor + 1
        testhalf <- drop(crossprod(stepPrev) < crossprod(step))
      }
      if (failedInv | (all(abs(step) < epsilon))) {
        break
      }
    }
  }

  ## ## If there was a failure in inverting the Fisher information then
  ## ## use the last invertible version of it, else recalculate everything
  ## if (failedInv) {
  ##   pars <- parsPrev
  ## }
  # OR
  ## Use always the most recent invertible version of the Fisher information
  pars <- parsPrev

  ## conduct single bias correction (if BC selected) else do not
  ## estimate the first order biases
  if (method == "BC") {
    bias <- biasFun(pars)$bias
    pars <- pars - bias
  }
  else {
    bias <- rep(NA_real_, pXL + q + pXS)
  }

  fit <- fitFun(pars)
  prob <- fit$prob
  scores <- gradFun(pars, fit = fit)
  infoInv <- infoFun(pars, fit = fit, inverse = TRUE)

  ## if method is BR then calcuate adjusted score functions else do
  ## not estimate their value
  if (method == "BR") {
    adjustedScores <- scores +
      biasFun(pars, fit = fit, vcov = infoInv)$adjustment
  }
  else {
     adjustedScores <- rep(NA_real_, pXL + q + pXS)
  }
  ## Just to cover the case where we revert back to the last value
  ## where the Fisher information was invertible
  failedInv <- failedInv | inherits(infoInv, "try-error")

  ## If maxit is not zero then reduce the number of iterations by 1,
  ## because we run till maxit + 1 and also pars <- parsPrev
  if (maxit > 0) {
    niter <- niter - 1
    if (history) {
        colnames(historyEstimates) <- c("Starting value",
        if (niter > 0) paste("Iteration", seq.int(length.out = niter)) else NULL)
    }
  }

  ## Extend the parameter vector to include aliased paraeters
  parsN <- structure(rep(NA, length(coefNames)), .Names = coefNames)
  parsN[!aliased] <- pars

  beta <- parsN[betaNames]
  alpha <- parsN[alphaNames]
  tau <- parsN[tauNames]

  ## Extend the inverse of the FIsher Information to include aliased parameters
  if (inherits(infoInv, "try-error")) {
    vcov <- infoInv
  }
  else {
    vcov <- structure(matrix(NA, length(coefNames), length(coefNames)),
                      .Dimnames = list(coefNames, coefNames))
    vcov[!aliased, !aliased] <- infoInv
  }

  inadmissible <- any((prob < 0) | (prob > 1))

  ## Return value
  rval <- list(beta = beta,
               alpha = alpha,
               tau = tau,
               bias = bias,
               scores = scores,
               adjustedScores = adjustedScores,
               vcov = vcov,
               call = match.call(),
               dataOriginal = if (missing(data)) NULL else data,
               dataReshaped = .dat,
               fitted.values = structure(c(prob), .Names = rownames(.dat)),
               deviance = -2 * sum(freqs[freqs != 0]*log(prob[freqs != 0])),
               convergence = (niter < maxit) & !failedInv,
               niter = niter,
               history = historyEstimates,
               boundaryEstimates = failedInv | inadmissible,
               inadmissible = inadmissible,
               method = method,
               eta = fit$eta,
               df.residual = sum(freqs) - pXL - q - pXS,
               edf = pXL + q + pXS,
               n = sum(freqs),
               nobs = sum(freqs), ## WHY?
               na.action = na.action,
               #xlevelsLocation = .getXlevels(TermsL, MLocation),
               model = if (model) mf else NULL,
               #xlevelsScale = if (pXS > 0) .getXlevels(TermsS, MScale) else NULL,
               contrasts = attr(XL, "contrasts"))
  class(rval) <- c("bpolr", "polr")
  rval
}

deviance.bpolr <- function(object) {
  object$deviance
}

evalFit.bpolr <- function(object, at) {
  method <- object$method
  update(object, start = at, method = if (method == "BC") "ML" else method,
         maxit = 0)
}

### Need print/summary methods that handle dispersion parameters
vcov.bpolr <- function(object, ...) {
  object$vcov
}


### Modify simulate.polr for handling dispersion parameters
simulate.bpolr <- function(object) {
  mf <- object$dataReshaped
  nlev <- nlevels(model.response(mf))
  weights <- model.weights(mf)
  N <- nrow(mf)
  NcovClass <- N/nlev
  probs <- matrix(object$fitted, nrow = nlev)
  totals <- tapply(weights, rep(1:NcovClass, each = nlev), sum)
  simFreqs <- vapply(1:length(totals),
                     function(i) rmultinom(1, totals[i], probs[,i]),
                     numeric(nlev))
  if (is.null(object$call$weights)) {
    mf["(weights)"] <- c(simFreqs)
  }
  else {
    mf["(weights)"] <- NULL
    weightsNam <- as.character(object$call$weights)
    mf[weightsNam] <- c(simFreqs)
  }
  mf
}


print.bpolr <- function (x, scoreDigits = 2, ...)
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    if (x$method == "BR") {
      cat("\nBias-reduced fit\n")
    }
    if (x$method == "ML") {
      cat("\nMaximum likelihood fit\n")
    }
    if (x$method == "BC") {
      cat("\nBias-corrected fit\n")
    }
    if (length(x$beta)) {
        cat("\nCoefficients:\n")
        print(x$beta, ...)
    }
    else {
        cat("\nNo coefficients\n")
    }
    cat("\nIntercepts:\n")
    print(x$alpha, ...)
    cat("\nScale parameters:\n")
    print(x$tau, ...)
    cat("\nResidual Deviance:", format(x$deviance, nsmall = 2L),
        "\n")
    cat("AIC:", format(x$deviance + 2 * x$edf, nsmall = 2L),
        "\n")
    if (nzchar(mess <- naprint(x$na.action)))
        cat("(", mess, ")\n", sep = "")
    if (x$boundaryEstimates) {
      cat("Warning: failed to invert the Fisher information. Iteration stopped prematurely\n")
    }
    if (x$inadmissible)
        cat("Warning: inadmissible parameter values\n")
    if (!x$convergence)
      cat("Warning: fitting algorithm did not converge\n")
    invisible(x)
}

coef.bpolr <- function(object,
                               what = c("all", "alpha", "beta", "tau")) {
  what <- match.arg(what)
  beta <- object$beta
  alpha <- object$alpha
  tau <- object$tau
  switch(what,
         all = c(beta, alpha, tau),
         alpha = alpha,
         beta = beta,
         tau = tau)
}

summary.bpolr <-
function (object, digits = max(3, .Options$digits - 3), correlation = FALSE,
    ...)
{
    cc <- c(object$beta, object$alpha, object$tau)
    pc <- length(object$beta)
    q <- length(object$alpha)
    pt <- length(object$tau)
    coef <- matrix(0, pc + q + pt, 3L, dimnames = list(names(cc),
        c("Value", "Std. Error", "t value")))
    coef[, 1L] <- cc
    vc <- vcov(object)
    if (inherits(vc, "try-error")) vc <- matrix(NA_real_, pc + q + pt,
                                                pc + q + pt)
    coef[, 2L] <- sd <- sqrt(diag(vc))
    coef[, 3L] <- coef[, 1L]/coef[, 2L]
    object$coefficients <- coef
    object$pc <- pc
    object$pt <- pt
    object$q <- q
    object$digits <- digits
    if (correlation)
        object$correlation <- (vc/sd)/rep(sd, rep(pc + q + pt, pc +
            q + pt))
    class(object) <- "summary.bpolr"
    object
}

print.summary.bpolr <- function (x, digits = x$digits, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl, control = NULL)
  }
  if (x$method == "BR") {
    cat("\nMaximum absolute adjusted score:\n")
    print(max(abs(x$adjustedScores)), ...)
    cat("\nIterations for bias-reduced fit: ", x$niter, "\n\n")
  }
  if (x$method == "ML") {
    cat("\nMaximum absolute score:\n")
    print(max(abs(x$scores)), ...)
    cat("\nIterations for maximum likelihood fit: ", x$niter, "\n\n")
  }
  if (x$method == "BC") {
    cat("\nBias-corrected fit\n\n")
  }
  if (x$pc > 0) {
    cat("\nCoefficients:\n")
    print(x$coefficients[seq_len(x$pc), , drop = FALSE], quote = FALSE,
          digits = digits, ...)
  }
  else {
    cat("\nNo coefficients\n")
  }
  cat("\nIntercepts:\n")
    print(x$coef[seq_len(x$q) + x$pc, , drop = FALSE], quote = FALSE,
          digits = digits, ...)
  if (x$pt > 0) {
    cat("\nScale parameters:\n")
    print(x$coef[seq_len(x$pt) + x$pc + x$q, , drop = FALSE], quote = FALSE,
          digits = digits, ...)
  }
  cat("\nResidual Deviance:", format(x$deviance, nsmall = 2L),
      "\n")
  cat("AIC:", format(x$deviance + 2 * x$edf, nsmall = 2L),
      "\n")
  if (nzchar(mess <- naprint(x$na.action)))
    cat("(", mess, ")\n", sep = "")
  if (!is.null(correl <- x$correlation)) {
    cat("\nCorrelation of Coefficients:\n")
    ll <- lower.tri(correl)
    correl[ll] <- format(round(correl[ll], digits))
    correl[!ll] <- ""
    print(correl[-1L, -ncol(correl)], quote = FALSE, ...)
  }
  invisible(x)
}
