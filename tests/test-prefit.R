require(fitdistrplus)

x <- rgamma(1e3, 5/2, 7/2)


fitdistrplus:::prefitmle(x, "dgamma", c(shape=3, scale=3), lower=-Inf, upper=Inf)

fitdistrplus:::prefitmle(x, "dgamma", c(shape=log(3), scale=log(3)), lower=0, upper=Inf)

fitdistrplus:::prefitmle(x, "dgamma", c(shape=log(3)), fix.arg=list(scale=7/2), lower=-Inf, upper=Inf)
