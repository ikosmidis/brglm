### Problems with subset
bpolr <- function(formula,
                  scale,
                  nominal,
                  data,
                  weights,
                  start = NULL,
                  subset,
                  na.action,
                  contrasts = NULL,
                  model = TRUE,
                  link = c("logit", "probit", "cloglog", "cauchit"),
                  method = c("ML", "BR", "BC"),
                  maxit = 100,
                  epsilon = 1e-10,
                  keepHistory = TRUE,
                  trace = FALSE,
                  slowIt = 1,
                  reparam = FALSE,
                  ...) {
  M  <- match.call(expand.dots = FALSE)
  if (missing(formula)) stop("A specification of location is absent")
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
  require(gnm)
  require(ordinal)
  method <- match.arg(method)
  if (is.matrix(eval.parent(M$data)))
    M$data <- as.data.frame(data)
  M$start <- M$method <- M$link <- M$model <- M$maxit <- M$epsilon <-
    M$keepHistory <- M$trace <- M$slowIt <- M$reparam <- M$... <- NULL
  M[[1L]] <- as.name("model.frame")
  MScale <- MLocation <- Mdata <- MNominal <- M
  MScale$formula <- MScale$nominal <-
      MLocation$scale <- MLocation$nominal <-
          Mdata$scale <- Mdata$formula <- Mdata$nominal <-
              MNominal$formula <- MNominal$scale <- NULL
  names(MScale)[names(MScale) == "scale"] <- "formula"
  names(MNominal)[names(MNominal) == "nominal"] <- "formula"
  respNam <- as.character(MLocation$formula[[2]])
  vars <- unique(c(all.vars(MLocation$formula),
                   all.vars(MScale$formula),
                   all.vars(MNominal$formula)))
  vars <- vars[vars!=respNam]
  ## Take care of intercept only models
  if (length(vars)) {
      ff <- as.formula(paste(respNam, "~", paste(vars, collapse = " + ")))
  }
  else
      ff <- as.formula(paste(respNam, "~ 1"))
  environment(ff) <- environment(formula)
  Mdata$formula <- ff
  .dat <- eval(Mdata)
  ## hack for intercept only models
  ## .dat$x <- c(1,1,1)
  termsMdata <- attr(.dat, "terms")
  nam <- names(.dat)
  if (!("(weights)"%in%nam)) {
    .dat[["(weights)"]] <- rep(1, nrow(.dat))
    nam <- c(nam, "(weights)")
  }
  MScale$weights <- MLocation$weights <- MNominal$weights <- as.name("(weights)")
  .dat <- expandCategorical(.dat, respNam, group = FALSE)
  .dat[["(weights)"]] <- c(.dat[["(weights)"]] * .dat[["count"]])
  inds <- match(c(nam[-1], nam[1]), names(.dat))
  .dat <- groupDat(.dat[, inds], respNam, "(weights)")
  .dat <- .dat[nam]
  attr(.dat, "terms") <- termsMdata
  MScale$data <- MLocation$data <- MNominal$data <- as.name(".dat")
  if (missing(scale)) MScale$formula <- ~ -1
  if (missing(nominal)) MNominal$formula <- ~ -1
  # A clm call for starting values later; used only if start is not specified
  Mclm <- M
  #  names(Mclm)[names(Mclm) == "location"] <- "formula"
  MLocation$subset <- MScale$subset <- MNominal$subset <- NULL
  #
  MLocation <- eval(MLocation)
  MScale <- eval(MScale)
  MNominal <- eval(MNominal)
  TermsL <- attr(MLocation, "terms")
  TermsS <- attr(MScale, "terms")
  TermsN <- attr(MNominal, "terms")
  # Useful integers
  nlev <- nlevels(model.response(MLocation))
  lev <- levels(model.response(MLocation))
  q <- nlev - 1
  N <- nrow(.dat)
  NcovClass <- N/nlev
  Nobs <- NcovClass*q
  rownames(.dat) <- paste(rep(1:NcovClass, each = nlev),
                          rep(1:nlev, times = NcovClass), sep = ".")
  ## Frequencies and totals
  freqs <- matrix(.dat[["(weights)"]], nrow = nlev)
  totals <- colSums(freqs)
  ## Set up model matrices
  X <- model.matrix(TermsL, MLocation)
  V <- model.matrix(TermsS, MScale)
  NN <- model.matrix(TermsN, MNominal)
  XX <- X[-seq(nlev, N, nlev), -1, drop = FALSE]
  .polr <- col(matrix(0, nrow(XX), q)) == rep(1:q, N/nlev)
  colnames(.polr) <- .polrNam <- paste(lev[-nlev], lev[-1L], sep = "|")
  ## Construct the nominal effects
  NN <- NN[-seq(nlev, N, nlev), -1, drop = FALSE]
  pN <- ncol(NN)
  nominalNames <- colnames(NN)
  NNew <- NULL
  if (pN > 0) {
      for (i in 1:pN) {
          Newmat <- NN[, i]*.polr
          colnames(Newmat) <- paste(.polrNam, nominalNames[i], sep = ".")
          NNew <- cbind(NNew, Newmat)
      }
  }
  XX <- cbind(XX, if (pN > 0) -NNew else NULL)
  Xlin <- cbind(-XX, .polr)
  # ADD NOMINAL SUPPORT
  VV <- V[-seq(nlev, N, nlev), -1, drop = FALSE]
  ## Set up offsets
  offsetL <- model.offset(MLocation)
  offsetS <- model.offset(MScale)
  if (is.null(offsetL)) offsetL <- rep(0, N)
  if (is.null(offsetS)) offsetS <- rep(0, N)
  ## offsetL will be substracted
  ## offsetV will be added
  offsetL <- offsetL[-seq(nlev, N, nlev)]
  offsetS <- offsetS[-seq(nlev, N, nlev)]
  # More useful integers
  pX <- ncol(XX)
  pV <- ncol(VV)
  inds <- sapply(1:NcovClass,
                 function(i) (1 + (i - 1)*q):(i*q))
  if (q > 1) {
    mitre <- diag(q) - cbind(0, rbind(diag(q - 1), 0))
  }
  else {
    dim(inds) <- c(1, NcovClass)
    mitre <- as.matrix(1)
  }
  ### Write etasExpr
  if (pV > 0) {
    colnames(VV) <- paste("scale", colnames(VV), sep = ".")
    etasExpr <- expression(matrix(c(drop(Xlin%*%pars[1:(pX + q)]) - offsetL)/
        exp(drop(VV%*%pars[pX + q + 1:pV] + offsetS)), nrow = q))
    Zexpr <- expression(cbind(Xlin/exp(drop(VV%*%pars[pX + q + 1:pV])),
        -VV*c(etas)))
  }
  else {
    etasExpr <- expression(matrix(c(drop(Xlin%*%pars[1:(pX + q)]) - offsetL),
        nrow = q))
    Zexpr <- expression(Xlin)
  }
  expInfo <- function(pars) {
    etas <- eval(etasExpr)
    cumprob <- pfun(etas)
    gs <- dfun(etas)
    prob <- apply(cumprob, 2, function(x) diff(c(0, x, 1)))
    FisherInfo <- 0
    Z <- eval(Zexpr)
    for (i in 1:NcovClass) {
      Zr <- Z[inds[,i], , drop = FALSE]
      invCovr <- 1/prob[nlev, i] + if (q == 1) 1/prob[1, i, drop = FALSE] else diag(1/prob[1:q, i])
      Dr <- mitre * gs[, i]
      Wr <- totals[i]*Dr%*%invCovr%*%t(Dr)
      FisherInfo <- FisherInfo + t(Zr)%*%Wr%*%Zr
    }
    FisherInfo
  }
  ders <- function(pars) {
    etas <- eval(etasExpr)
    cumprob <- pfun(etas)
    gs <- dfun(etas)
    prob <- apply(cumprob, 2, function(x) diff(c(0, x, 1)))
    FisherInfo <- 0
    Z <- eval(Zexpr)
    for (i in 1:NcovClass) {
      Zr <- Z[inds[,i], , drop = FALSE]
      invCovr <- 1/prob[nlev, i] + if (q == 1) 1/prob[1, i, drop = FALSE] else diag(1/prob[1:q, i])
      Dr <- mitre * gs[, i]
      Wr <- totals[i]*Dr%*%invCovr%*%t(Dr)
      FisherInfo <- FisherInfo + t(Zr)%*%Wr%*%Zr
    }
    list(scores = -colSums(c(gs*apply(freqs/prob, 2, diff))*Z),
         expInfo = FisherInfo)
  }
  # Need adjustment for dispersion pars
  adjustment <- function(pars) {
    etas <- eval(etasExpr)
    cumprob <- pfun(etas)
    gs <- dfun(etas)
    prob <- apply(cumprob, 2, function(x) diff(c(0, x, 1)))
    Z <- eval(Zexpr)
    invExpInfo <- solve(expInfo(pars))
    cc <- matrix(0, nrow(prob), ncol(prob))
    for (i in 1:NcovClass) {
      etar <- etas[, i]
      Zr <- Z[inds[,i], , drop = FALSE]
      Ar <- Zr%*%invExpInfo%*%t(Zr)
      cr <- 0.5*totals[i]*diag(Ar)*ddfun(etar)
      cc[, i] <- c(cr[1], diff(cr), -cr[q])
    }
    -colSums(c(gs*apply(cc/prob, 2, diff))*Z)
  }
  adjustmentDisp <- function(pars) {
    etas <- eval(etasExpr)
    cumprob <- pfun(etas)
    gs <- dfun(etas)
    prob <- apply(cumprob, 2, function(x) diff(c(0, x, 1)))
    Z <- eval(Zexpr)
    invExpInfo <- solve(expInfo(pars))
    invFtt <- invExpInfo[pX + q + 1:pV, pX + q + 1:pV]
    invFat <- invExpInfo[pX + 1:q, pX + q + 1:pV]
    invFbt <- invExpInfo[1:pX, pX + q + 1:pV]
    taus <- pars[q + pX + 1:pV]
    adj <- adj1 <- adj2 <- 0
    cc <- matrix(0, nrow(prob), ncol(prob))
    for (i in 1:NcovClass) {
      etar <- etas[, i]
      invCovr <- 1/prob[nlev, i] + if (q == 1) 1/prob[1, i, drop = FALSE] else diag(1/prob[1:q, i])
      Dr <- mitre * gs[, i]
      Wr <- totals[i]*Dr%*%invCovr%*%t(Dr)
      Zr <- Z[inds[, i], , drop = FALSE]
      vr <- VV[inds[1, i], ]
      xr <- XX[inds[1, i], ]
      Ar <- Zr%*%invExpInfo%*%t(Zr)
      cr <- 0.5*totals[i]*diag(Ar)*ddfun(etar)
      cc[, i] <- c(cr[1], diff(cr), -cr[q])
      adj <- adj + 0.5*(vr%*%invFtt%*%vr) * c(Wr%*%etar)*Zr -
        exp(-sum(taus*vr))*c(Wr%*%invFat%*%vr -
                             if (pX > 0) xr%*%invFbt%*%vr * rowSums(Wr)
                             else 0)*Zr
    }
    colSums(adj) - colSums(c(gs*apply(cc/prob, 2, diff))*Z)
  }
  BC <- FALSE
  if (method == "BC") {
    maxit <- slowIt <- 1
    BC <- TRUE
    method <- "BR"
    tempCall <- match.call()
    tempCall$method <- "ML"
    tempObj <- eval(tempCall)
    pars <- c(tempObj$beta, tempObj$alpha, tempObj$tau)
    ### Dropped dependence on clm for the MLE here
    ## browser()
    ## www <- .dat[["(weights)"]]
    ## Mclm[[1L]] <- as.name("clm")
    ## Mclm$link <- as.name("link")
    ## Mclm$weights <-as.name("www")
    ## Mclm$maxIter <- 100
    ## datS <- .dat
    ## datS$www <- www
    ## datS$iii <- www > 0
    ## Mclm$data <- as.name("datS")
    ## Mclm$subset <- as.name("iii")
    ## ## options(warn = -1)
    ## clmObject <- eval(Mclm)
    ## ## options(warn = 0)
    ## ## zeta is what is called tau here
    ## pars <- c(clmObject$beta, clmObject$alpha[-c(1:q)],
    ##           clmObject$alpha[1:q], clmObject$zeta)
  }
  else {
      if (is.null(start)) {
          ## Starting values
          ## Add something small to zero weights to avoid the
          ## possibility of boundary estimates
          www <- .dat[["(weights)"]]
          www[www == 0] <- 0.5*(pX + q + pV)/sum(www)
          Mclm[[1L]] <- as.name("clm")
          Mclm$link <- as.name("link")
          Mclm$weights <-as.name("www")
          Mclm$maxIter <- 2
          datS <- .dat
          datS$www <- www
          Mclm$data <- as.name("datS")
          options(warn = -1)
          clmObject <- eval(Mclm)
          options(warn = 0)
          ## zeta is what is called tau here
          pars <- c(clmObject$beta, clmObject$alpha[-c(1:q)],
                    clmObject$alpha[1:q], clmObject$zeta)
      }
      else {
          if (missing(scale)) {
              if (length(start)!=(pX + q))
                  stop("'start' is not of the correct length", call. = FALSE)
          }
          else {
              if (length(start)!=(pX + pV + q))
                  stop("'start' is not of the correct length", call. = FALSE)
          }
          pars <- start
      }
  }
  niter <- 0
  historyEsts <- pars
  infEsts <- FALSE
  if (missing(scale)) {
    bias <- function(pars) {
      -solve(expInfo(pars))%*%adjustment(pars)
    }
  }
  else {
    bias <- function(pars) {
      -solve(expInfo(pars))%*%adjustmentDisp(pars)
    }
  }
  if (method == "BR") {
    if (missing(scale)) {
      gradExpr <- expression(curDers$scores + adjustment(pars))
    }
    else {
      gradExpr <- expression(curDers$scores + adjustmentDisp(pars))
    }
  }
  else {
    gradExpr <- expression(curDers$scores)
  }
  ### Pars is in beta, alpha, tau
  ### Pars1 is in beta, theta, tau
  JacAlphaTheta <- function(pars1) {
    theta <- pars1[pX + 1:q]
    etheta <- exp(theta)
    mat <- matrix(0, q, q)
    mat[, 1L] <- rep(1, q)
    if (q > 1)
      for (i in 2L:q) mat[i:q, i] <- etheta[i]
    J <- diag(pX + q + pV)
    J[pX + 1:q, pX + 1:q] <- mat
    J
  }
  dersTheta <- function(pars1) {
    theta <- pars1[pX + 1:q]
    alpha <- cumsum(c(theta[1], exp(theta[-1])))
    parsOR <- pars1
    parsOR[pX + 1:q] <- alpha
    dersOR <- ders(parsOR)
    A <- JacAlphaTheta(pars1)
    scoresTheta <- dersOR$scores%*%A
    infoTheta <- t(A)%*%dersOR$expInfo%*%A
    list(scores = structure(c(dersOR$scores%*%A),
           names = names(dersOR$scores)),
         expInfo = structure(t(A)%*%dersOR$expInfo%*%A,
           dimnames = dimnames(dersOR$expInfo)))
  }
  biasTheta <- function(pars1) {
    theta <- pars1[pX + 1:q]
    alpha <- cumsum(c(theta[1], exp(theta[-1])))
    parsOR <- pars1
    parsOR[pX + 1:q] <- alpha
    invFisherInfo <- solve(expInfo(parsOR))
    biasPars <- bias(parsOR)
    biasAlpha <- biasPars[pX + 1:q]
    invFaa <- invFisherInfo[pX + 1:q, pX + 1:q]
    if (q == 1) {
      biasPars[pX + 1] <- biasAlpha[1]
    }
    else {
      ttt <- sapply(2:q, function(i)
                    sum(c(-1, 1, 1, -1)*c(invFaa[(i-1):i, (i-1):i])))
      biasPars[pX + 1:q] <- c(biasAlpha[1],
                              diff(biasAlpha)/diff(alpha) + 0.5*ttt/diff(alpha)^2)
    }
    biasPars
  }
  if (reparam) {
    pars1 <- pars
    alpha <- pars[pX + 1:q]
    theta <- c(alpha[1], log(diff(alpha)))
    pars1[pX + 1:q] <- theta
    while (niter < maxit) {
      curDers <- dersTheta(pars1)
      curInfo <- curDers$expInfo
      invCurInfo <- try(solve(curInfo), silent = TRUE)
      if (inherits(invCurInfo, "try-error")) {
        infEsts <- TRUE
        break
      }
      niter <- niter + 1
      grad <- curDers$scores - if (method == "BR") curInfo%*%biasTheta(pars1) else 0
      if (trace) {
        cat("===========================\n")
        cat("Inverse condition number:", 1/kappa(curInfo, exact = TRUE), "\n")
        cat("Scores:\n", format(grad, scientific = TRUE), "\n")
        cat("===========================\n")
        browser()
      }
      pars1 <- pars1 + slowIt*invCurInfo%*%grad
      if (keepHistory) historyEsts <- cbind(historyEsts, pars1)
      if (all(abs(grad) < epsilon)) break
    }
    pars <- pars1
    theta <- pars1[pX + 1:q]
    alpha <- cumsum(c(theta[1], exp(theta[-1])))
    pars[pX + 1:q] <- alpha
  }
  else {
  ### A quasi Fisher scoring iteration
    while (niter < maxit) {
      curDers <- ders(pars)
      curInfo <- curDers$expInfo
      invCurInfo <- try(solve(curInfo), silent = TRUE)
      if (inherits(invCurInfo, "try-error")) {
        infEsts <- TRUE
        break
      }
      niter <- niter + 1
      grad <- eval(gradExpr)
      if (trace) {
        cat("===========================\n")
        cat("Inverse condition number:", 1/kappa(curInfo, exact = TRUE), "\n")
        cat("Scores:\n", format(grad, scientific = TRUE), "\n")
        cat("===========================\n")
        browser()
      }
      pars <- pars + slowIt*invCurInfo%*%grad
      if (keepHistory) historyEsts <- cbind(historyEsts, pars)
      if (all(abs(grad) < epsilon)) break
    }
  }
  parNames <- rownames(pars)
  pars <- c(pars)
  names(pars) <- parNames
  cumprob <- pfun(eval(etasExpr))
  prob <- apply(cumprob, 2, function(x) diff(c(0, x, 1)))
  curDers <- ders(pars)
  object <- NULL
  if (pX!=0) {
    object$beta <- pars[1:pX]
    object$alpha <- pars[pX + 1:q]
  }
  else {
    object$alpha <- pars[1:q]
  }
  if (!missing(scale)) {
    object$tau <- pars[pX + q + 1:pV]
  }
  object$expInfo <- curDers$expInfo
  if (any(is.na(curDers$expInfo))) {
    infEsts <- TRUE
  }
  dimnames(object$expInfo) <- list(parNames, parNames)
  object$scores <- curDers$scores
  if (infEsts) {
    warning("Either some of the maximum likelihood estimates are on the boundary  \n of the parameter space or better starting values need to be supplied", call. = FALSE)
    object$adjustedScores <- rep(NA, pX + q + pV)
#    object$adjustedScores <- grad
    if (method == "ML") object$firstOrderBias <- rep(NA, pX + q + pV)
  }
  else {
    if (method == "BR") {
      if (reparam) {
        curDers <- dersTheta(pars1)
        object$adjustedScores <- curDers$scores - curDers$expInfo%*%biasTheta(pars1)
      }
      else {
        if (missing(scale)) {
          object$adjustedScores <- curDers$scores + adjustment(pars)
        }
        else {
          object$adjustedScores <- curDers$scores + adjustmentDisp(pars)
        }
      }
    }
    else {
      # because first order biases do not make sense when method == "BR"
      object$firstOrderBias <- bias(pars)
      pars1 <- pars
      alpha <- pars[pX + 1:q]
      theta <- c(alpha[1], log(diff(alpha)))
      pars1[pX + 1:q] <- theta
      object$firstOrderBiasReparam <- biasTheta(pars1)
    }
  }
  object$pars <- pars
  object$call <- match.call()
  object$dataReshaped <- .dat
  object$fitted.values <- c(prob)
  names(object$fitted.values) <- rownames(.dat)
  object$deviance <- -2*sum(freqs*log(prob))
  object$convergence <- as.numeric(((niter >= maxit)*!BC) || infEsts)
  object$niter <- niter
  if (!is.null(historyEsts) & !infEsts) {
    colnames(historyEsts) <- c("Starting", paste("Iteration", 1:niter))
  }
  object$history <- historyEsts
  object$infiniteEstimates <- infEsts
  object$method <- if (BC) "BC" else method
  object$eta <- eval(etasExpr)
  object$df.residual <- sum(freqs) - pX - q - pV
  object$edf <- pX + q + pV
  object$n <- sum(freqs)
  object$nobs <- sum(freqs)
  object$na.action <- attr(MLocation, "na.action")
  object$xlevelsLocation <- .getXlevels(TermsL, MLocation)
  object$reparam <- reparam
  if (model)
    object$model <- .dat
  if (pV > 0)
    object$xlevelsScale <- .getXlevels(TermsS, MScale)
  object$contrasts <- attr(X, "contrasts")
  class(object) <- c("bpolr", "polr")
  object
}

model.frame.bpolr <- function(formula, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots),
                    0)]
  if (length(nargs) || is.null(formula$model)) {
    M  <- formula$call
    if (is.matrix(eval.parent(M$data)))
      M$data <- as.data.frame(data)
    M$start <- M$method <- M$link <- M$model <- M$maxit <- M$epsilon <-
      M$keepHistory <- M$trace <- M$slowIt <- M$... <- NULL
    M[[1L]] <- as.name("model.frame")
    MScale <- MLocation <- Mdata <- M
    MScale$formula <- MLocation$scale <- Mdata$scale <- Mdata$formula <-NULL
    names(MScale)[names(MScale) == "scale"] <- "formula"
    respNam <- as.character(MLocation$formula[[2]])
    vars <- unique(c(all.vars(MLocation$formula), all.vars(MScale$formula)))
    vars <- vars[vars!=respNam]
    ff <- as.formula(paste(respNam, "~", paste(vars, collapse = " + ")))
    respNam <- as.character(MLocation$formula[[2]])
    Mdata$formula <- ff
    .dat <- eval(Mdata)
    termsMdata <- attr(.dat, "terms")
    nam <- names(.dat)
    if (!("(weights)"%in%nam)) {
      .dat[["(weights)"]] <- rep(1, nrow(.dat))
      nam <- c(nam, "(weights)")
    }
    MScale$weights <- MLocation$weights <- as.name("(weights)")
    .dat <- expandCategorical(.dat, respNam, group = FALSE)
    .dat[["(weights)"]] <- c(.dat[["(weights)"]] * .dat[["count"]])
    inds <- match(c(nam[-1], nam[1]), names(.dat))
    .dat <- groupDat(.dat[, inds], respNam, "(weights)")
    .dat <- .dat[nam]
    attr(.dat, "terms") <- termsMdata
    .dat
  }
  else formula$model
}

### Need print/summary methods that handle dispersion parameters
vcov.bpolr <- function(object, ...) {
  solve(object$expInfo)
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
    if (x$convergence > 0)
        cat("Warning: did not converge as iteration limit reached\n")
    invisible(x)
}

coef.bpolr <- function(object, ...) {
    c(object$beta, object$alpha, object$ta)
}

coefficients.bpolr <- function(object, ...) {
    coef(object, ...)
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
    print(x$coef[(x$pc + 1L):(x$pc + x$q), , drop = FALSE], quote = FALSE,
          digits = digits, ...)
  cat("\nScale parameters:\n")
  print(x$coef[(x$pc + x$q + 1L):(x$pc + x$q + x$pt), , drop = FALSE], quote = FALSE,
        digits = digits, ...)
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
