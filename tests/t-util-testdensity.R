library(fitdistrplus)


testdpqfun <- fitdistrplus:::testdpqfun

##### existence ##### 
#a list of TRUE
testdpqfun("exp")
#a list of error messages
testdpqfun("exp2")

##### void vector ##### 
dexp2 <- function(x, rate)
  ifelse(length(x)==0, stop("zero input"), dexp(x,rate))
dexp3 <- function(x, rate)
  ifelse(length(x)==0, NA, dexp(x,rate))
#TRUE
testdpqfun("exp", "d", c(rate=1))
#error message
testdpqfun("exp2", "d", c(rate=1))
#error message
testdpqfun("exp3", "d", c(rate=1))

##### inconsistent value ##### 
pexp2 <- function(q, rate)
{
  res <- pexp(q, rate)
  if(any(is.nan(res)))
    stop("NaN values")
  res
}
pexp3 <- function(q, rate)
{
  res <- pexp(q, rate)
  if(any(is.infinite(q)))
    stop("Inf values")
  res
}

#TRUE
testdpqfun("exp", "p", c(rate=1))
#error message
testdpqfun("exp2", "p", c(rate=1))
#error message
testdpqfun("exp3", "p", c(rate=1))

##### missing value ##### 
qexp2 <- function(p, rate)
{
  res <- qexp(p, rate)
  if(any(is.na(res)))
    stop("NA values")
  res
}
qexp3 <- function(p, rate)
{
  res <- qexp(p, rate)
  res[!is.na(res)]
}

#TRUE
testdpqfun("exp", "q", c(rate=1))
#error message
testdpqfun("exp2", "q", c(rate=1))
#error message
testdpqfun("exp3", "q", c(rate=1))

##### inconsistent parameter ##### 
dnorm2 <- function(x, mean, sd)
{
  if(sd < 0)
    stop("negative param")
  else
    dnorm(x,mean,sd)
}
#TRUE
testdpqfun("norm", "d", c(mean=1, sd=1))
#error message
testdpqfun("norm2", "d", c(mean=1, sd=1))

##### inconsistent name ##### 
dnorm2 <- function(x, mean=0, sd=1, ...)
    dnorm(x,mean,sd)

#TRUE
testdpqfun("norm", "d", c(mean=1, sd=1))
#error message
testdpqfun("norm2", "d", c(mean=1, sd=1))
