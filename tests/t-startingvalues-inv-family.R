
require(actuar)
require(fitdistrplus)
n <- 1e4


#### Inverse exponential ####
set.seed(123)
x <- rinvexp(n, 3)

if(FALSE)
{
  fitdistrplus:::startarg_invtransgamma_family(x, "invexp")
  fitdistrplus:::startarg_transgamma_family(1/x, "exp")
}

try(cdfcomp(fitdist(x, "invexp", lower=1e-6), xlogscale = TRUE, do.points = FALSE))


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


try(cdfcomp(fitdist(x, "invtrgamma", lower=1e-6, upper=100, start=list("shape1"=50, "shape2"=1, "scale"=1/2)), 
            xlogscale = TRUE, do.points = FALSE))

x <- rinvparalogis(n, 2, 10)


try(cdfcomp(fitdist(x, "invparalogis", lower=1e-6, start=list("shape"=1e-5, "scale"=1/2)), 
            xlogscale = TRUE, do.points = FALSE))

