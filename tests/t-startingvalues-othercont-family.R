library(fitdistrplus)
library(actuar)

set.seed(1234)

n <- 1e3

x <- rlgamma(n, 1, 1)
fitdistrplus:::startarg_othercontinuous_actuar_family(x, "lgamma")
fitdist(x, "lgamma")


x <- rgumbel(n, 1, 1)
fitdistrplus:::startarg_othercontinuous_actuar_family(x, "gumbel")
fitdist(x, "gumbel")


x <- rinvgauss(n, pi, 2*pi)
fitdistrplus:::startarg_othercontinuous_actuar_family(x, "invgauss")
fitdist(x, "invgauss")


x <- rgenbeta(n, 2, 2, 2, scale=2)
fitdistrplus:::startarg_othercontinuous_actuar_family(x, "genbeta")
fitdist(x, "genbeta", lower=0)
