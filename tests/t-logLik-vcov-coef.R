library(fitdistrplus)

# (1) basic fit of a gamma distribution by maximum likelihood estimation
#
data(groundbeef)
serving <- groundbeef$serving
fitg <- fitdist(serving, "gamma")

logLik(fitg)
vcov(fitg)
coef(fitg)

fitg <- fitdist(serving, "gamma", method="mme")

logLik(fitg)
vcov(fitg) 
coef(fitg)

# (2) Fit of a lognormal distribution to bacterial contamination data
#
data(smokedfish)
fitsf  <-  fitdistcens(smokedfish,"lnorm")
logLik(fitsf)
vcov(fitsf) 
coef(fitsf)

