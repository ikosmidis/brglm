\name{separation.detection}
\alias{separation.detection}
\title{Separation Identification.}
\description{
  Provides a tool for identifying whether or not separation has occurred.
}
\usage{
separation.detection(fit, nsteps = 30)
}
\arguments{
  \item{fit}{the result of a \code{\link{glm}} call.}
  \item{nsteps}{Starting from \code{maxit = 1}, the GLM is refitted for
  \code{maxit = 2}, \code{maxit = 3}, \ldots, \code{maxit = nsteps}. Default
  value is 30.}
}
\details{
  Identifies separated cases for binomial-response GLMs, by refitting
  the model. At each iteration the maximum number of allowed IWLS
  iterations is fixed starting from 1 to \code{nsteps} (by setting
  \code{control = glm.control(maxit = j)}, where \code{j} takes values  1,
  \ldots, nsteps in \code{\link{glm}}). For each value of \code{maxit},
  the estimated asymptotic standard errors are divided to the
  corresponding ones resulted for
  \code{control = glm.control(maxit = 1)}. Based on the results in Lesaffre
  & Albert (1989), if the sequence of ratios in any column of the
  resultant matrix diverges, then separation occurs and the maximum
  likelihood estimate for the corresponding parameter has value minus or
  plus infinity.
}
\value{
  A matrix of dimension \code{nsteps} by \code{length(coef(fit))}, that
  contains the ratios of the estimated asymptotic standard errors.
}
\references{
Lesaffre, E. and Albert, A. (1989). Partial separation in logistic
discrimination. \emph{J. R. Statist. Soc. \bold{B}}, \bold{51}, 109--116.
}

\author{Ioannis Kosmidis, \email{i.kosmidis@ucl.ac.uk}}

\examples{
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
}
\keyword{models}
\keyword{utilities}