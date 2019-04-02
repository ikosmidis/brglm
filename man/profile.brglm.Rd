\name{profile.brglm}
\alias{profile.brglm}
\alias{print.profile.brglm}
\title{Calculate profiles for objects of class 'brglm'.}
\description{
  Creates \code{"profile.brglm"} objects to be used for the calculation of
  confidence intervals and for plotting.
}
\usage{
\method{profile}{brglm}(fitted, gridsize = 10, stdn = 5,
        stepsize = 0.5, level = 0.95,
        which = 1:length(coef(fitted)), verbose = TRUE,
        zero.bound = 1e-08, scale = FALSE, ...)
}
\arguments{
  \item{fitted}{an object of class \code{"brglm"}.}
  \item{gridsize}{as in \code{\link[profileModel]{profileModel}}.}
  \item{stdn}{as in \code{\link[profileModel]{profileModel}}.}
  \item{stepsize}{as in \code{\link[profileModel]{profileModel}}.}
  \item{level}{\code{qchisq(level,1)} indicates the range that the
    profiles must cover.}
  \item{which}{as in \code{\link[profileModel]{profileModel}}.}
  \item{verbose}{as in \code{\link[profileModel]{profileModel}}.}
  \item{zero.bound}{as in \code{\link[profileModel]{profileModel}}.}
  \item{scale}{as in \code{\link[profileModel]{profileModel}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  \code{profile.brglm} calculates the profiles of the appropriate
  objectives to be used for the construction of confidence intervals for
  the bias-reduced estimates (see \code{\link{confint.brglm}} for the
  objectives that are profiled).
}
\value{
  An object of class \code{"profile.glm"} with attribute \dQuote{level}
  corresponding to the argument \code{level}. The object supports the
  methods \code{\link{print}}, \code{\link{plot}}, \code{\link{pairs}}
  and \code{\link{confint}} and it  is a list of the components:
  \item{profilesML}{a \code{"profileModel"} object containing
    the profiles of the ordinary deviance for the maximum likelihood fit
    corresponding to \code{fitted}.}
  \item{profilesBR}{\code{NULL} if \code{method = "glm.fit"} in
  \code{\link{brglm}}. If \code{method = "brglm.fit"} and \code{pl =
  TRUE}, \code{profilesBR} is a \code{"profileModel"} object containing
  the profiles of the penalized deviance for the parameters of
  \code{fitted}. If \code{method = "brglm.fit"} and \code{pl = FALSE}
  \code{profilesBR} is a \code{"profileModel"} object containing the
  profiles of the modified score statistic (see
  \code{\link{profileObjectives}}) for the parameters of \code{fitted}.}
}
\note{
  Objects of class \code{"profile.brglm"} support the methods:
  \itemize{
    \item{\code{print}}{which prints the \code{"level"} attribute of the
      object, as well as the supported methods.}
    \item{\code{confint}}{see \code{\link{confint.brglm}}.}
    \item{\code{plot}}{see \code{\link{plot.profile.brglm}}.}
    \item{\code{pairs}}{see \code{\link{plot.profile.brglm}}.}
  }
}
\author{Ioannis Kosmidis, \email{i.kosmidis@ucl.ac.uk}}
\seealso{\code{\link[profileModel]{profileModel}}, \code{\link{profile.brglm}}.}
\examples{
# see example in 'confint.brglm'.
}
\keyword{models}