require("fitdistrplus")
require("actuar")

set.seed(1234)

n <- 1e3
x <- rztpois(n, 2)

fitdistrplus:::startarg_discrete_actuar_family(x, "ztpois")
fitdist(x, "ztpois")

fitdistrplus:::startarg_discrete_actuar_family(x, "ztnbinom")
fitdist(x, "ztnbinom")

fitdistrplus:::startarg_discrete_actuar_family(x, "ztgeom")
fitdist(x, "ztgeom")

x <- rztbinom(n, 30, 1/2)

fitdistrplus:::startarg_discrete_actuar_family(x, "ztbinom")
fitdist(x, "ztbinom", lower=0)
