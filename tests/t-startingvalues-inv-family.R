
require(actuar)
require(fitdistrplus)
n <- 1e4


#### Inverse exponential ####

x <- rinvexp(n, 3)

if(FALSE)
{
fitdistrplus:::startarg_invtransgamma_family(x, "invexp")
fitdistrplus:::startarg_transgamma_family(1/x, "exp")
}

cdfcomp(fitdist(x, "invexp"), xlogscale = TRUE, do.points = FALSE)


#### Inverse transformed gamma ####

x <- rinvtrgamma(n, 3, 3, 10)

if(FALSE)
{
  fitdistrplus:::startarg_invtransgamma_family(x, "invtrgamma")
  fitdistrplus:::startarg_transgamma_family(1/x, "trgamma")
  
  cutshapeparam(list("shape0"=1e3, "shape33"=1e-10, "theta"=22))
  cutshapeparam(list("shape0"=1e3, "shape33"=1e3, "theta"=22))
  
  cutshapeparam(fitdistrplus:::startarg_invtransgamma_family(x, "invtrgamma"))
}

cdfcomp(fitdist(x, "invtrgamma", lower=0), xlogscale = TRUE, do.points = FALSE)
cdfcomp(fitdist(x, "invtrgamma", lower=0, start=list("shape1"=100, "shape2"=1, "scale"=1/2)), 
        xlogscale = TRUE, do.points = FALSE)
curve(pinvtrgamma(x, 10, 1, 1))
curve(pinvtrgamma(x, 25, 1, 1/2))
curve(pinvtrgamma(x, 100, 1, 1/2))

x <- rinvparalogis(n, 2, 10)

cdfcomp(fitdist(x, "invparalogis"), xlogscale = TRUE, do.points = FALSE)
cdfcomp(fitdist(x, "invparalogis", lower=0, start=list("shape"=1e-5, "scale"=1/2),
                control=list(trace=1, REPORT=1)), 
        xlogscale = TRUE, do.points = FALSE)

