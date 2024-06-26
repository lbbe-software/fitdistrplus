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


#weird examples
x <- rep(1, n)

fitdistrplus:::startarg_transgamma_family(x, "gamma")

#previous code
n <- length(x)
m <- mean(x)
v <- (n - 1)/n*var(x)
list(shape=m^2/v, rate=m/v)


#normal -> weibull
x <- abs(rnorm(n = 20, mean = 150, sd = 10))
m <- mean(log(x))
v <- var(log(x))
shape <- 1.2/sqrt(v)
scale <- exp(m + 0.572/shape)
s1 <- list(shape=shape, scale=scale)
cdfcomp(fitdist(x, "weibull", start=s1))

