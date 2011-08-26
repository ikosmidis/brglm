\name{confint.brglm}
\alias{confint.brglm}
\alias{confint.profile.brglm}
\title{Computes confidence intervals of parameters for bias-reduced estimation}
\description{
  Computes confidence intervals for one or more parameters when
  estimation is performed using \code{\link{brglm}}. The resultant confidence
  intervals are based on manipulation of the profiles of the deviance,
  the penalized deviance and the modified score statistic (see
  \code{\link{profileObjectives}}).
}
\usage{
\method{confint}{brglm}(object, parm = 1:length(coef(object)), level = 0.95,
        verbose = TRUE, endpoint.tolerance = 0.001,
        max.zoom = 100, zero.bound = 1e-08, stepsize = 0.5,
        stdn = 5, gridsize = 10, scale = FALSE, method = "smooth",
        ci.method = "union", n.interpolations = 100, ...)

\method{confint}{profile.brglm}(object, parm, level = 0.95, method = "smooth",
        ci.method = "union", endpoint.tolerance = 0.001,
        max.zoom = 100, n.interpolations = 100, verbose = TRUE,
        ...)
}
\arguments{
  \item{object}{an object of class \code{"brglm"} or \code{"profile.brglm"}.}
  \item{parm}{either a numeric vector of indices or a character vector
    of names, specifying the parameters for which confidence intervals
    are to be estimated. The default is all parameters in the fitted
    model. When \code{object} is of class \code{"profile.brglm"},
    \code{parm} is not used and confidence intervals are returned for
    all the parameters for which profiling took place.}
  \item{level}{the confidence level required. The default is 0.95. When
    \code{object} is of class \code{"profile.brglm"}, \code{level} is not
    used and the level attribute of \code{object} is used instead.}
  \item{verbose}{logical. If \code{TRUE} (default) progress indicators
  are printed during the progress of calculating the confidence
  intervals.}
  \item{endpoint.tolerance}{as in \code{\link[profileModel]{confintModel}}.}
  \item{max.zoom}{as in \code{\link[profileModel]{confintModel}}.}
  \item{zero.bound}{as in \code{\link[profileModel]{confintModel}}.}
  \item{stepsize}{as in \code{\link[profileModel]{confintModel}}.}
  \item{stdn}{as in \code{\link[profileModel]{confintModel}}.}
  \item{gridsize}{as in \code{\link[profileModel]{confintModel}}.}
  \item{scale}{as in \code{\link[profileModel]{confintModel}}.}
  \item{method}{as in \code{\link[profileModel]{confintModel}}.}
  \item{ci.method}{The method to be used for the construction of
    confidence intervals. It can take values \code{"union"} (default) and
    \code{"mean"} (see Details).}
  \item{n.interpolations}{as in \code{\link[profileModel]{confintModel}}.}
  \item{\dots}{further arguments to or from other methods.}
}
\details{
  In the case of logistic regression Heinze & Schemper (2002) and Bull
  et. al. (2007) suggest the use of confidence intervals based on the
  profiles of the penalized likelihood, when estimation is performed
  using maximum penalized likelihood.

  Kosmidis (2007) illustrated that because of the shape of the
  penalized likelihood, confidence intervals based on the penalized
  likelihood could exhibit low or even zero coverage for hypothesis
  testing on large parameter values and also misbehave illustrating severe
  oscillation (see Brown et. al., 2001). Kosmidis (2007)  suggested an
  alternative confidence interval that is based on the union of the
  confidence intervals  resulted by profiling the ordinary deviance for
  the maximum likelihood fit and by profiling the penalized deviance for
  the maximum penalized fit.  Such confidence intervals,  despite of
  being slightly conservative,  illustrate less oscillation and avoid
  the loss of coverage. Another possibility is to use the mean of the
  corresponding endpoints instead of \dQuote{union}. Yet unpublished
  simulation studies suggest that such confidence intervals are not as
  conservative as the \dQuote{union} based intervals but illustrate more
  oscillation, which however is not as severe as in the case of the
  penalized likelihood based ones.
  
  The properties of the \dQuote{union} and \dQuote{mean} confidence
  intervals extend  to all the links that are supported by
  \code{\link{brglm}}, when  estimation is performed using maximum
  penalized likelihood. 

  In the case of estimation using modified scores and for models other
  than logistic, where there is not an objective that is maximized, the
  profiles of the penalized likelihood for the construction of the
  \dQuote{union} and \dQuote{mean} confidence intervals can be replaced
  by the profiles of modified score statistic (see
  \code{\link{profileObjectives}}).

  The \code{confint} method for \code{brglm} and \code{profile.brglm}
  objects implements the \dQuote{union} and \dQuote{mean} confidence
  intervals. The method is chosen through the \code{ci.method} argument.
}
\value{
  A matrix with columns the endpoints of the confidence intervals for
  the specified (or profiled) parameters.
}
\references{
  Brown, L. D., Cai, T. T. and DasGupta, A. (2001).
  Interval estimation for a binomial proportion (with discussion).
 \emph{Statistical Science} \bold{16}, 101--117.
  
  Bull, S. B., Lewinger, J. B. and Lee, S. S. F. (2007). Confidence
  intervals for multinomial logistic regression in sparse data.
  \emph{Statistics in Medicine} \bold{26}, 903--918.

  Heinze, G. and Schemper, M. (2002). A solution to the problem of
  separation in logistic regression.
  \emph{Statistics in Medicine} \bold{21}, 2409--2419.

  Kosmidis, I. (2007). Bias reduction in exponential family nonlinear
  models. \emph{PhD Thesis}, Department of Statistics, University of
  Warwick.
}
\author{Ioannis Kosmidis, \email{ioannis@stats.ucl.ac.uk}}

\seealso{\code{\link[profileModel]{confintModel}}, \code{\link[profileModel]{profileModel}},
  \code{\link{profile.brglm}}. }
\examples{
## Begin Example 1
\dontrun{
#library(MASS)
#data(bacteria)
#contrasts(bacteria$trt) <- structure(contr.sdif(3),
#          dimnames = list(NULL, c("drug", "encourage")))
# fixed effects analyses
#m.glm.logit <- brglm(y ~ trt * week, family = binomial,
#                     data = bacteria, method = "glm.fit")
#m.brglm.logit <- brglm(y ~ trt * week, family = binomial,
#                       data = bacteria, method = "brglm.fit")
#p.glm.logit <- profile(m.glm.logit)
#p.brglm.logit <- profile(m.brglm.logit)
# 
#plot(p.glm.logit)
#plot(p.brglm.logit)
# confidence intervals for the glm fit based on the profiles of the
# ordinary deviance
#confint(p.glm.logit)
# confidence intervals for the brglm fit
#confint(p.brglm.logit, ci.method = "union")
#confint(p.brglm.logit, ci.method = "mean")
# A cloglog link
#m.brglm.cloglog <- update(m.brglm.logit, family = binomial(cloglog))
#p.brglm.cloglog <- profile(m.brglm.cloglog)
#plot(p.brglm.cloglog)
#confint(m.brglm.cloglog, ci.method = "union")
#confint(m.brglm.cloglog, ci.method = "mean")
## End example
## Begin Example 2
y <- c(1, 1, 0, 0)
totals <- c(2, 2, 2, 2)
x1 <- c(1, 0, 1, 0)
x2 <- c(1, 1, 0, 0)
m1 <- brglm(y/totals ~ x1 + x2, weights = totals,
            family = binomial(cloglog))
p.m1 <- profile(m1)
confint(p.m1, method="zoom")
}
}
\keyword{models}
\keyword{htest}
