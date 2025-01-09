
require(actuar)
require(fitdistrplus)
n <- 1e4

x <- rinvexp(n, 3)

if(FALSE)
{
fitdistrplus:::startarg_invtransgamma_family(x, "invexp")
fitdistrplus:::startarg_transgamma_family(1/x, "exp")
}

cdfcomp(fitdist(x, "invexp"), xlogscale = TRUE, do.points = FALSE)


x <- rinvtrgamma(n, 3, 3, 10)

if(FALSE)
{
  fitdistrplus:::startarg_invtransgamma_family(x, "invtrgamma")
  fitdistrplus:::startarg_transgamma_family(1/x, "exp")
}

cdfcomp(fitdist(x, "invtrgamma", lower=0), xlogscale = TRUE, do.points = FALSE)



x <- rinvparalogis(n, 2, 10)

cdfcomp(fitdist(x, "invparalogis"), xlogscale = TRUE, do.points = FALSE)
