require(fitdistrplus)
data("fremale")

# fremale test

fremale.cens <- Surv2fitdistcens(fremale$AgeIn, fremale$AgeOut, fremale$Death)


f1 <- fitdistcens(fremale.cens, "norm")
f2 <- fitdistcens(fremale.cens, "logis")
f3 <- fitdistcens(fremale.cens, "cauchy")

cdfcompcens(list(f1, f2, f3))
sapply(list(f1, f2, f3), logLik)

