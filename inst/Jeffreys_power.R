## Ioannis Kosmidis, 11 Sep 2019

library("brglm")

## Set up custom additive modifications to the responses and totals;
## see ?brglm::modifications and the examples there for details

## .const below is the multiplier of the log-determinant of the Fisher information:
## .const = 1/2 does bias reduction
## brglm will retrieve the value of .const from the global environment here
## (not neat but it works with sufficient care!)
br.custom.family <- function(p) {
     list(ar = .const * p/p, at = 2 * .const * p/p)
}

## Set up a custom link-glm object (essentially just copying logit
## here as we only care about modifying the penalty for logit links)
mylogit <- make.link("logit")
mylogit$name <- "mylogit"

data("lizards")

## The reduced-bias fit is
.const <- 1/2
brglm(cbind(grahami, opalinus) ~ height + diameter +
          light + time, family = binomial(mylogit), data=lizards)

## which is the same as what brglm does by default for logistic regression
brglm(cbind(grahami, opalinus) ~ height + diameter +
          light + time, family = binomial(logit), data=lizards)


## Stronger penalization (e.g. 5/2) can be achieved by 
.const <- 5/2
brglm(cbind(grahami, opalinus) ~ height + diameter +
          light + time, family = binomial(mylogit), data=lizards)

