source("brglm.fit.R")
source("brglm.control.R")
source("confint.brglm.R")
source("myfamily.R")
source("plot.profile.brglm.R")
source("profileObjectives.R")
source("separation.detection.R")


### An example from
### King, Gary, James E. Alt, Nancy Elizabeth Burns and Michael Laver
### (1990).  "A Unified Model of Cabinet Dissolution in Parliamentary
### Democracies", _American Journal of Political Science_, vol. 34,
### no. 3, pp. 846-870.

library(Zelig)
data(coalition)
fm1 <- glm(duration ~ fract + numst2, family = Gamma, data = coalition)
clotting <- data.frame(u = c(5,10,15,20,30,40,60,80,100),
                       lot1 = c(118,58,42,35,27,25,21,19,18),
                       lot2 = c(69,35,26,21,18,16,13,12,12))

## The bias-reduced fit. The bias-reduced estimate of the dispersion
## parameter is also calculated.
fm1br <- update(m1, method = "brglm.fit")
## The bias-corrected fit. The bias-corrected estimate of the dispersion
## parameter is also calculated. The correction argument is from brglm.control
fm1bc <- update(m1, method = "brglm.fit", correction = TRUE)

## Now another bias-reduced fit. Bias is reduced for 1/dispersion
## instead of dispersion
fm1brInv <- glm(duration ~ fract + numst2, family = Gamma, data = coalition,
                method = "brglm.fit", dispTrans = "inverse")
## Now another bias-corrected fit. Bias is corrected for 1/dispersion
## instead of dispersion
fm1bcInv <- glm(duration ~ fract + numst2, family = Gamma, data = coalition,
                method = "brglm.fit", dispTrans = "inverse", correction = TRUE)




### The lizards example from the previous version of brglm
data(lizards)
## Fit the model using maximum likelihood
lizardsML <- glm(cbind(grahami, opalinus) ~ height + diameter +
                 light + time, family = binomial(logit), data=lizards,
                 method = "glm.fit")
## Now the bias-reduced fit:
lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter +
                 light + time, family = binomial(logit), data=lizards,
                 method = "brglm.fit")
lizardsML
≈lizardsBR
## Other links
update(lizardsBR, family = binomial(probit))
update(lizardsBR, family = binomial(cloglog))
update(lizardsBR, family = binomial(cauchit))
