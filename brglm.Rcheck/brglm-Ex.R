pkgname <- "brglm"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('brglm')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("brglm")
### * brglm

flush(stderr()); flush(stdout())

### Name: brglm
### Title: Bias reduction in Binomial-response GLMs
### Aliases: brglm brglm.fit print.brglm summary.brglm print.summary.brglm
### Keywords: models regression iteration

### ** Examples

## Begin Example
data(lizards)
# Fit the GLM using maximum likelihood
lizards.glm <- brglm(cbind(grahami, opalinus) ~ height + diameter +
                  light + time, family = binomial(logit), data=lizards,
                  method = "glm.fit")
# Now the bias-reduced fit:
lizards.brglm <- brglm(cbind(grahami, opalinus) ~ height + diameter +
                  light + time, family = binomial(logit), data=lizards,
                  method = "brglm.fit")
lizards.glm
lizards.brglm
# Other links
update(lizards.brglm, family = binomial(probit))
update(lizards.brglm, family = binomial(cloglog))
update(lizards.brglm, family = binomial(cauchit))
# Using penalized maximum likelihood
update(lizards.brglm, family = binomial(probit), pl = TRUE)
update(lizards.brglm, family = binomial(cloglog), pl = TRUE)
update(lizards.brglm, family = binomial(cauchit), pl = TRUE)



cleanEx()
nameEx("confint.brglm")
### * confint.brglm

flush(stderr()); flush(stdout())

### Name: confint.brglm
### Title: Computes confidence intervals of parameters for bias-reduced
###   estimation
### Aliases: confint.brglm confint.profile.brglm
### Keywords: models htest

### ** Examples

## Begin Example 1
## Not run: 
##D library(MASS)
##D data(bacteria)
##D contrasts(bacteria$trt) <- structure(contr.sdif(3),
##D           dimnames = list(NULL, c("drug", "encourage")))
##D # fixed effects analyses
##D m.glm.logit <- brglm(y ~ trt * week, family = binomial,
##D                      data = bacteria, method = "glm.fit")
##D m.brglm.logit <- brglm(y ~ trt * week, family = binomial,
##D                        data = bacteria, method = "brglm.fit")
##D p.glm.logit <- profile(m.glm.logit)
##D p.brglm.logit <- profile(m.brglm.logit)
##D # 
##D plot(p.glm.logit)
##D plot(p.brglm.logit)
##D # confidence intervals for the glm fit based on the profiles of the
##D # ordinary deviance
##D confint(p.glm.logit)
##D # confidence intervals for the brglm fit
##D confint(p.brglm.logit, ci.method = "union")
##D confint(p.brglm.logit, ci.method = "mean")
##D # A cloglog link
##D m.brglm.cloglog <- update(m.brglm.logit, family = binomial(cloglog))
##D p.brglm.cloglog <- profile(m.brglm.cloglog)
##D plot(p.brglm.cloglog)
##D confint(m.brglm.cloglog, ci.method = "union")
##D confint(m.brglm.cloglog, ci.method = "mean")
##D ## End example
## End(Not run)
## Begin Example 2
y <- c(1, 1, 0, 0)
totals <- c(2, 2, 2, 2)
x1 <- c(1, 0, 1, 0)
x2 <- c(1, 1, 0, 0)
m1 <- brglm(y/totals ~ x1 + x2, weights = totals,
            family = binomial(cloglog))
p.m1 <- profile(m1)
confint(p.m1, method="zoom")



cleanEx()
nameEx("lizards")
### * lizards

flush(stderr()); flush(stdout())

### Name: lizards
### Title: Habitat Preferences of Lizards
### Aliases: lizards
### Keywords: datasets

### ** Examples

data(lizards)
glm(cbind(grahami, opalinus) ~ height + diameter + light + time,
    family = binomial, data=lizards)
brglm(cbind(grahami, opalinus) ~ height + diameter + light + time,
    family = binomial, data=lizards)



cleanEx()
nameEx("modifications")
### * modifications

flush(stderr()); flush(stdout())

### Name: modifications
### Title: Additive Modifications to the Binomial Responses and Totals for
###   Use within 'brglm.fit'
### Aliases: modifications checkModifications
### Keywords: models regression

### ** Examples

## Begin Example 1
## logistic exposure model, following the Example in ?family. See,
## Shaffer, T.  2004. Auk 121(2): 526-540.
# Definition of the link function
logexp <- function(days = 1) {
  linkfun <- function(mu) qlogis(mu^(1/days))
  linkinv <- function(eta) plogis(eta)^days
  mu.eta <- function(eta) days * plogis(eta)^(days-1) *
        .Call("logit_mu_eta", eta, PACKAGE = "stats")
  valideta <- function(eta) TRUE
  link <- paste("logexp(", days, ")", sep="")
  structure(list(linkfun = linkfun, linkinv = linkinv,
    mu.eta = mu.eta, valideta = valideta, name = link),
    class = "link-glm")
}
# Here d(p) = days * p * ( 1 - p^(1/days) )
#      d'(p) = (days - (days+1) * p^(1/days)) * d(p)
#      w(p) = days^2 * p * (1-p^(1/days))^2 / (1-p)
# Initial modifications, as given from the general expressions above:
br.custom.family <- function(p) {
  etas <- binomial(logexp(.days))$linkfun(p)
  # the link function argument `.days' will be detected by lexical
  # scoping. So, make sure that the link-function inputted arguments
  # have unusual names, like `.days' and that
  # the link function enters `brglm' as
  # `family=binomial(logexp(.days))'. 
  list(ar=0.5*(1-p)-0.5*(1-p)*exp(etas)/.days,
       at=0*p/p) # so that to fix the length of at
}
.days <-3
# `.days' could be a vector as well but then it should have the same
# length as the number of observations (`length(.days)' should be
# equal to `length(p)'). In this case, `checkModifications' should
# have argument `Length=length(.days)'.
#
# Check:
## Not run: checkModifications(br.custom.family)
# OOOPS error message... the condition is not satisfied
#   
# After some trivial algebra using the two allowed operations, we
# get new modifications:
br.custom.family <- function(p) {
  etas <- binomial(logexp(.days))$linkfun(p)
  list(ar=0.5*p/p, # so that to fix the length of ar
       at=0.5+exp(etas)*(1-p)/(2*p*.days))
}
# Check:
checkModifications(br.custom.family)
# It is OK.
# Now,
modifications(binomial(logexp(.days)))
# works.
# Notice that for `.days <- 1', `logexp(.days)' is the `logit' link
# model and `a_r=0.5', `a_t=1'.
# In action:
library(MASS)
example(birthwt)
m.glm <- glm(formula = low ~ ., family = binomial, data = bwt)
.days <- bwt$age
m.glm.logexp <- update(m.glm,family=binomial(logexp(.days)))
m.brglm.logexp <- brglm(formula = low ~ ., family =
binomial(logexp(.days)), data = bwt)
# The fit for the `logexp' link via maximum likelihood
m.glm.logexp
# and the fit for the `logexp' link via modified scores
m.brglm.logexp
## End Example
## Begin Example 2
## Another possible use of brglm.fit:
## Deviating from bias reducing modified-scores:
## Add 1/2 to the response of a probit model.
y <- c(1,2,3,4)
totals <- c(5,5,5,5)
x1 <- c(1,0,1,0)
x2 <- c(1,1,0,0)
my.probit <- make.link("probit")
my.probit$name <- "my.probit"
br.custom.family <- function(p) {
   h <- pmax(get("hats",parent.frame()),.Machine$double.eps)
   list(ar=0.5/h,at=1/h)
}
m1 <- brglm(y/totals~x1+x2,weights=totals,family=binomial(my.probit))
m2 <- glm((y+0.5)/(totals+1)~x1+x2,weights=totals+1,family=binomial(probit))
# m1 and m2 should be the same.    



cleanEx()
nameEx("plot.profile.brglm")
### * plot.profile.brglm

flush(stderr()); flush(stdout())

### Name: plot.profile.brglm
### Title: Plot methods for 'profile.brglm' objects
### Aliases: plot.profile.brglm pairs.profile.brglm
### Keywords: dplot hplot

### ** Examples

# see example in 'confint.brglm'.



cleanEx()
nameEx("profile.brglm")
### * profile.brglm

flush(stderr()); flush(stdout())

### Name: profile.brglm
### Title: Calculate profiles for objects of class 'brglm'.
### Aliases: profile.brglm print.profile.brglm
### Keywords: models

### ** Examples

# see example in 'confint.brglm'.



cleanEx()
nameEx("separation.detection")
### * separation.detection

flush(stderr()); flush(stdout())

### Name: separation.detection
### Title: Separation Identification.
### Aliases: separation.detection
### Keywords: models utilities

### ** Examples

## Begin Example
y <- c(1,1,0,0)
totals <- c(2,2,2,2)
x1 <- c(1,0,1,0)
x2 <- c(1,1,0,0)
m1 <- glm(y/totals ~ x1 + x2, weights = totals, family = binomial())
# No warning from glm...
m1
# However estimates for (Intercept) and x2 are unusually large in
# absolute value... Investigate further:
#
separation.detection(m1,nsteps=30)
# Note that the values in the column for (Intercept) and x2 diverge,
# while for x1 converged. Hence, separation has occurred and the
# maximum lieklihood estimate for (Intercept) is minus infinity and
# for x2 is plus infinity. The signs for infinity are taken from the
# signs of (Intercept) and x1 in coef(m1).
## End Example



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
