require(fitdistrplus)

#case where the density returns a Not-an-Numeric value.
detectbound("gamma", c(shape=3, scale=3), 1:10)
detectbound("logis", c(location=3, scale=3), 1:10)
detectbound("geom", c(prob=1/2), 1:10)

#test rate-scale arg
detectbound("gamma", c(shape=1, scale=3), 1:10)
detectbound("gamma", c(shape=1, rate=1/3), 1:10)


#case where the density returns a Not-an-Numeric value and one parameter is fixed.
detectbound("gamma", c(shape=3), 1:10, fix.arg=c(scale=3))


#case where the density returns an error rather than a Not-an-Numeric value.
dgeom2 <- function(x, prob, log=FALSE)
{
  stopifnot(prob >= 0 && prob <= 1)
  dgeom(x, prob, log)
}
detectbound("geom2", c(prob=1/2), 1:10)

#case where the density returns a Not-an-Numeric value for actuar package
require(actuar)

detectbound("burr", c(shape1=3, shape2=3, rate=1), 1:10)
detectbound("llogis", c(shape=3, rate=1), 1:10)
