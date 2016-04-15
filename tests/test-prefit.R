require(fitdistrplus)

x <- rgamma(1e3, 5/2, 7/2)

fitdistrplus:::prefitmle(x, "gamma", "mle", c(shape=3, scale=3), lower=-Inf, upper=Inf, silent=TRUE, control=list(trace=1, REPORT=1))
fitdistrplus:::prefitmle(x, "gamma", "mle", c(shape=1, scale=1), lower=-Inf, upper=Inf, silent=TRUE)

fitdistrplus:::prefitmle(x, "gamma", "mle", c(shape=3), fix.arg=list(scale=7/2), lower=-Inf, upper=Inf, silent=TRUE)


x <- rgeom(1e3, 1/7)

fitdistrplus:::prefitmle(x, "geom", "mle", c(prob=1/2), lower=-Inf, upper=Inf, silent=TRUE)

tbx <- table(x)
fitdistrplus:::prefitmle(as.numeric(names(tbx)), "geom", "mle", c(prob=1/2), lower=-Inf, upper=Inf, silent=TRUE, weights=tbx)

1/7
