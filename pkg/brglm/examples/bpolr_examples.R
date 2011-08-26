rm(list = ls())
source("bpolr.R")
library(MASS)
library(ordinal)
library(gnm)

## The standard example in ?polr
data(housing)
house.plr <- polr(formula = Sat ~ Infl + Type + Cont,
                  data = housing, weights = Freq,
                  Hess = TRUE)
## Let's randomize the order of the data and remove some observations
houseDatNew <- housing[sample(1:72, 72, replace = FALSE), ]
houseDatNew <- houseDatNew[-c(1,2,7, 11, 33, 69), ]
houseDatNew$Freq <- round(houseDatNew$Freq/6)
## Fit a model with main effects Infl, Type and Cont and scale effects Infl.
fm1ML <- bpolr(formula = Sat ~ Infl + Type + Cont,
                scale = ~ Infl,
                data = houseDatNew,
                weights = Freq,
                link = "logit",
                method = "ML")
## Check against clm to make sure we are doing the right thing
fm2ML <- clm(formula = Sat ~ Infl + Type + Cont,
              scale = ~ Infl,
              data = houseDatNew,
              weights = Freq,
              link = "logit",
              subset = Freq > 0)
fm1BR <- update(fm1ML, method = "BR")
fm1BC <- update(fm1ML, method = "BC")
summary(fm1ML)
## Call:
## bpolr(formula = Sat ~ Infl + Type + Cont, scale = ~Infl, data = houseDatNew,
##     weights = Freq, link = "logit", method = "ML")
##
## Scores:
##       InflMedium         InflHigh    TypeApartment       TypeAtrium
##     8.902219e-13    -4.826556e-13    -2.132590e-11     2.830902e-12
##      TypeTerrace         ContHigh       Low|Medium      Medium|High
##    -6.954215e-12    -9.684462e-12    -1.796854e-12    -3.906077e-13
## scale.InflMedium   scale.InflHigh
##     2.273163e-11     6.912800e-12
##
## Iterations for maximum likelihood fit:  23
##
##
## Coefficients:
##                 Value Std. Error t value
## InflMedium     0.4831     0.2793   1.730
## InflHigh       1.9045     0.6595   2.888
## TypeApartment -0.7902     0.3569  -2.214
## TypeAtrium    -0.4953     0.4433  -1.117
## TypeTerrace   -1.4365     0.4680  -3.069
## ContHigh       0.2889     0.2726   1.060
##
## Intercepts:
##               Value Std. Error t value
## Low|Medium  -0.8545     0.3947  -2.165
## Medium|High  0.4362     0.3685   1.184
##
## Scale parameters:
##                    Value Std. Error t value
## scale.InflMedium 0.08765     0.2286  0.3834
## scale.InflHigh   0.44990     0.3404  1.3216
##
## Residual Deviance: 544.9587
## AIC: 564.9587
summary(fm1BR)
## Call:
## bpolr(formula = Sat ~ Infl + Type + Cont, scale = ~Infl, data = houseDatNew,
##     weights = Freq, link = "logit", method = "BR")
##
## Adjusted scores:
##       InflMedium         InflHigh    TypeApartment       TypeAtrium
##    -1.418421e-12     3.276268e-13     6.731171e-12    -8.473222e-13
##      TypeTerrace         ContHigh       Low|Medium      Medium|High
##     1.360967e-12     3.288703e-12     2.922773e-12    -2.445821e-13
## scale.InflMedium   scale.InflHigh
##    -3.261835e-11    -9.652390e-12
##
## Iterations for bias-reduced fit:  23
##
##
## Coefficients:
##                 Value Std. Error t value
## InflMedium     0.4658     0.2772   1.680
## InflHigh       1.7513     0.6177   2.835
## TypeApartment -0.7548     0.3524  -2.142
## TypeAtrium    -0.4739     0.4386  -1.080
## TypeTerrace   -1.3774     0.4623  -2.980
## ContHigh       0.2814     0.2699   1.043
##
## Intercepts:
##               Value Std. Error t value
## Low|Medium  -0.8162     0.3906  -2.090
## Medium|High  0.4322     0.3654   1.183
##
## Scale parameters:
##                    Value Std. Error t value
## scale.InflMedium 0.07597     0.2316  0.3281
## scale.InflHigh   0.41910     0.3381  1.2397
##
## Residual Deviance: 545.0762
## AIC: 565.0762
summary(fm1BC)
## Call:
## bpolr(formula = Sat ~ Infl + Type + Cont, scale = ~Infl, data = houseDatNew,
##     weights = Freq, link = "logit", method = "BC")
##
## Bias-corrected fit
##
##
## Coefficients:
##                 Value Std. Error t value
## InflMedium     0.4665     0.2784   1.676
## InflHigh       1.7473     0.6175   2.830
## TypeApartment -0.7654     0.3542  -2.161
## TypeAtrium    -0.4798     0.4405  -1.089
## TypeTerrace   -1.3909     0.4644  -2.995
## ContHigh       0.2798     0.2709   1.033
##
## Intercepts:
##               Value Std. Error t value
## Low|Medium  -0.8282     0.3923  -2.111
## Medium|High  0.4254     0.3666   1.161
##
## Scale parameters:
##                    Value Std. Error t value
## scale.InflMedium 0.08411     0.2315  0.3633
## scale.InflHigh   0.42393     0.3375  1.2560
##
## Residual Deviance: 545.087
## AIC: 565.087
