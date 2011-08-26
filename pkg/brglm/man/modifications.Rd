\name{modifications}

\alias{modifications}

\alias{checkModifications}

\title{Additive Modifications to the Binomial Responses and Totals for
Use within `brglm.fit'}

\description{
  Get, test and set the functions that calculate the additive
  modifications to the responses and totals in binomial-response GLMs,
  for the application of bias-reduction either via modified scores or
  via maximum penalized likelihood  (where penalization is by Jeffreys
  invariant prior).
}

\usage{
modifications(family, pl = FALSE)
}

\arguments{
  \item{family}{a family object of the form \code{binomial(link = "link")}, where
  \code{"link"} can be one of \code{"logit"}, \code{"probit"}, \code{"cloglog"}
  and \code{"cauchit"}. The usual ways of giving the family name are
  supported (see \code{\link{family}}).}
  \item{pl}{logical determining whether the function returned corresponds
  to modifications for the penalized maximum likelihood approach or for
  the modified-scores approach to bias-reduction. Default value is
  \code{FALSE}.}
}

\details{
  The function returned from \code{modifications} accepts the argument \code{p}
  which are the binomial probabilities and returns a list with
  components \code{ar} and \code{at}, which are the link-dependent parts
  of the additive modifications to the actual responses and totals,
  respectively.

  Since the resultant function is used in \code{\link{brglm.fit}}, for
  efficiency reasons no check is made for \code{p >= 0 | p <= 1}, for
  \code{length(at) == length(p)} and for \code{length(ap) == length(p)}.
}

\section{Construction of custom pseudo-data representations}{
  If \eqn{y^*}{y*} are the pseudo-responses (pseudo-counts) and
  \eqn{m^*}{m*} are the pseudo-totals then we call the pair \eqn{(y^*,
  m^*)}{(y*, m*)} a pseudo-data representation. Both the modified-scores approach and the maximum penalized likelihood have a common property:

there exists \eqn{(y^*, m^*)}{(y*, m*)} such that if we replace the actual data
  \eqn{(y, m)} with \eqn{(y^*, m^*)}{(y*, m*)} in the expression for the
  ordinary scores (first derivatives of the likelihood) of a
  binomial-response GLM, then we end-up either with the modified-scores
  or with the derivatives of the penalized likelihood (see Kosmidis,
  2007, Chapter 5).

  Let \eqn{\mu} be the mean of the binomial response \eqn{y}
  (i.e. \eqn{\mu=mp}{\mu = m p}, where \eqn{p} is the binomial probability
  corresponding to the count \eqn{y}). Also, let \eqn{d} and \eqn{d'}
  denote the first and the second derivatives, respectively, of
  \eqn{\mu}{\mu} with respect to the linear predictor \eqn{\eta}{\eta} of the
  model. All the above are viewed as functions of \eqn{p}. The
  pseudo-data representations have the generic form
  
  \tabular{ll}{  
    pseudo-response : \tab \eqn{y^*=y + h a_r(p)}{y* = y + h a_r(p)} \cr
    pseudo-totals : \tab \eqn{m^*=m + h a_t(p)}{m* = m + h a_t(p)}, \cr
  }

  where \eqn{h} is the leverage corresponding to \eqn{y}. The general
  expressions for \eqn{a_r(p)} ("r" for "response") and \eqn{a_t(p)}
  ("t" for "totals") are:

  \emph{modified-scores approach}
  \tabular{l}{
   \eqn{a_r(p) = d'(p)/(2w(p))} \cr
   \eqn{a_t(p) = 0}, \cr
  }

  \emph{maximum penalized likelihood approach}
  \tabular{l}{
   \eqn{a_r(p) = d'(p)/w(p) + p - 0.5} \cr
   \eqn{a_t(p) = 0}. \cr
  }
  
  For supplying \eqn{(y^*, m^*)}{(y*, m*)} in \code{\link{glm.fit}} (as is
  done by \code{\link{brglm.fit}}), an essential requirement for the
  pseudo-data representation is that it should mimic the behaviour of the
  original responses and totals, i.e. \eqn{0 \le y^* \le m^*}{0 \le y*
  \le m*}. Since \eqn{h \in [0, 1]}, the requirement translates to
  \eqn{0 \le a_r(p) \le  a_t(p)} for every \eqn{p \in (0, 1)}. However,
  the above definitions of \eqn{a_r(p)} and \eqn{a_t(p)} do not
  necessarily respect this requirement.

  On the other hand, the pair \eqn{(a_r(p), a_t(p))} is not unique in
  the sense that for a given link function and once the link-specific
  structure of the pair has been extrapolated, there is a class of
  equivalent pairs that can be resulted following only the following two
  rules:

  \itemize{
    \item add and subtract the same quantity from either \eqn{a_r(p)}
    or \eqn{a_t(p)}.
    \item if a quantity is to be moved from \eqn{a_r(p)} to \eqn{a_t(p)} it
    first has to be divided by \eqn{-p}.
  }
  
  For example, in the case of penalized maximum likelihood, the pairs
  \eqn{(d'(p)/w(p) + p - 0.5 , 0)} and \eqn{(d'(p)/w(p) + p , 0.5/p)} are
  equivalent, in the sense that if the corresponding pseudo-data
  representations are substituted in the ordinary scores both return the
  same expression.
 
  So, in order to construct a pseudo-data representation that
  corresponds to a user-specified link function and has the property
  \eqn{0 \le a_r(p) \le  a_t(p)} for every \eqn{p \in (0, 1)}, one merely
  has to pursue a simple algebraic calculation on the initial pair
  \eqn{(a_r(p), a_t(p))} using only the two aforementioned rules until
  an appropriate pair is resulted. There is always a pair!

  Once the pair has been found the following steps should be followed.
  \enumerate{
    \item For a user-specified link function the user has to write a
          modification function with name "br.custom.family" or
          "pml.custom.family" for \code{pl=FALSE} or \code{pl=TRUE},
          respectively. The function should take as argument the
          probabilities \code{p} and return a list  of two vectors with
          same length as \code{p} and  with names
          \code{c("ar", "at")}. The result corresponds to the pair
          \eqn{(a_r(p), a_t(p))}.
	  \item Check if the custom-made modifications function is
	  appropriate. This can be done via the function
          \code{\link{checkModifications}} which has arguments
          \code{fun} (the function to be tested) and \code{Length} with
          default value \code{Length=100}. \code{Length} is to be used
          when the user-specified link function takes as argument a
          vector of values (e.g. the \code{logexp} link in
          \code{?family}). Then the value of \code{Length} should be the
          length of that vector.
    \item Put the function in the search patch so that
          \code{modifications} can find it.
    \item \code{\link{brglm}} can now be used with the custom family as
          \code{\link{glm}} would be used.
  }
}
\note{
  The user could also deviate from modified-scores and maximum penalized
  likelihood and experiment with implemented (or not) links, e.g. \code{probit},
  constructing his own pseudo-data representations of the aforementioned
  general form. This could be done by changing the link name, e.g. by

  \code{probitt <-  make.link(probit) ;
    probitt$name <- "probitt"}

  and then setting a custom \code{br.custom.family} that does
  not necessarily depend on the \code{probit} link. Then, \code{brglm}
  could be used with \code{pl=FALSE}.

  A further generalization would be to completely remove the hat value
  \eqn{h} in the generic expression of the pseudo-data representation
  and have general additive modifications that depend on \eqn{p}. To do
  this divide both \code{ar} and \code{at} by
  \code{pmax(get("hats",parent.frame()),.Machine\$double.eps)} within the
  custom modification function (see also Examples).
}
\author{Ioannis Kosmidis, \email{ioannis@stats.ucl.ac.uk}}
\references{
  Kosmidis, I. (2007). Bias reduction in exponential family nonlinear
  models. \emph{PhD Thesis}, Department of Statistics, University of
  Warwick.
}
\seealso{\code{\link{brglm}}, \code{\link{brglm.fit}}}
\examples{
\dontrun{
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
\dontrun{checkModifications(br.custom.family)}
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
}
}
\keyword{models}
\keyword{regression}
