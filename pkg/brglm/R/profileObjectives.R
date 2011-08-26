adjustedScoreStatistic <- function (fm, X, dispersion = 1)
{
  # check for zero weights and zero varmu?
    y <- fm$y
    wt <- fm$prior.weights
    LP <- fm$linear.predictor
    family <- fm$family
    mu <- family$linkinv(LP)
    dmu <- family$mu.eta(LP)
    ddmu <- family$dmu.deta(LP)
    varmu <- family$variance(mu)
    w <- wt * dmu^2/varmu
    W.X <- sqrt(w) * X
    XWXinv <- chol2inv(chol(crossprod(W.X)))
    hats <- diag(X %*% XWXinv %*% t(w * X))
    scoresComp <- wt * dmu/varmu * (y - mu) * X
    adjExpComp <- 0.5 * hats * ddmu/dmu * X
    scores <- colSums(scoresComp)
    adjExp <- colSums(adjExpComp)
    adjscores <- scores/dispersion + adjExp
    dispersion * t(adjscores) %*% XWXinv %*% adjscores
}


### Needs work...
penalizedDeviance <- function (fm, X, dispersion = 1)
{
    Y <- fm$y
    LP <- fm$linear.predictor
    fam <- fm$family
    wt <- fm$prior.weights
    mu <- fm$fitted.values
    we <- fm$weights
    W.X <- sqrt(we) * X
    (sum(fam$dev.resid(Y, mu, wt)) - log(det(crossprod(W.X))))/dispersion
}
