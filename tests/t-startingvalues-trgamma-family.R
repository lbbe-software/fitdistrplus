library(fitdistrplus)
library(actuar)

set.seed(1234)

n <- 1e3
x <- rtrgamma(n, 2, 2, scale=2)

fitdistrplus:::startarg_transgamma_family(x, "trgamma")

fitdistrplus:::startarg_transgamma_family(x, "gamma")

fitdistrplus:::startarg_transgamma_family(x, "weibull")

fitdistrplus:::startarg_transgamma_family(x, "exp")

fitdist(x, "trgamma")
fitdist(x, "gamma")
fitdist(x, "weibull")
fitdist(x, "exp")


data(groundbeef)
fitdist(groundbeef$serving, "gamma", fix.arg = list(rate = 0.5))
