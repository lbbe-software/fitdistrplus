library(fitdistrplus)
library(actuar)

set.seed(1234)

n <- 1e3
x <- rzmpois(n, 2, 1/2)

fitdistrplus:::startarg_discrete_actuar_family(x, "zmpois")
fitdist(x, "zmpois")

x <- rzmnbinom(n, 2, 1/2, 1/3)
fitdistrplus:::startarg_discrete_actuar_family(x, "zmnbinom")
fitdist(x, "zmnbinom")

x <- rzmgeom(n, 1/2, 1/3)
fitdistrplus:::startarg_discrete_actuar_family(x, "zmgeom")
fitdist(x, "zmgeom")


x <- rzmbinom(n, 30, 1/2, 1/3)

fitdistrplus:::startarg_discrete_actuar_family(x, "zmbinom")
fitdist(x, "zmbinom", lower=0)


x <- rpoisinvgauss(n, mean = 10, dispersion = 2)

fitdistrplus:::startarg_discrete_actuar_family(x, "poisinvgauss")
fitdist(x, "poisinvgauss", lower=0)

