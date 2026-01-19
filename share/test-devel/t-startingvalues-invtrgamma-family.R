require("fitdistrplus")
require("actuar")

set.seed(1234)

n <- 1e3
x <- rinvtrgamma(n, 2, 2, scale=2)

fitdistrplus:::startarg_invtransgamma_family(x, "invtrgamma")

fitdistrplus:::startarg_invtransgamma_family(x, "invgamma")

fitdistrplus:::startarg_invtransgamma_family(x, "invweibull")

fitdistrplus:::startarg_invtransgamma_family(x, "invexp")

fitdist(x, "invtrgamma")
fitdist(x, "invgamma")
fitdist(x, "invweibull")
fitdist(x, "invexp")

