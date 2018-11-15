library(fitdistrplus)


testdpqfun <- fitdistrplus:::testdpqfun



##### first argument ##### 
#a data.frame of TRUE and ""
testdpqfun("exp", start=c(rate=1))
#a data.frame with error messages
dEXP <- function(y, rate) dexp(x, rate)
pEXP <- function(y, rate) pexp(x, rate)
qEXP <- function(y, rate) qexp(x, rate)
testdpqfun("EXP", start=c(rate=1))


##### existence ##### 
#a data.frame of TRUE and ""
testdpqfun("exp", start=c(rate=1))
#a data.frame with error messages
testdpqfun("exp2", start=c(rate=1))

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
    dnorm(x, mean, sd)

dnorm3 <- dnorm2
pnorm3 <- pnorm
qnorm3 <- qnorm

#TRUE
testdpqfun("norm", "d", c(mean=1, sd=1))
#error message
testdpqfun("norm2", "d", c(mean=1, sd=1))


#a data.frame with error messages
testdpqfun("norm", c("d", "p", "q"), c(mean=1, sd=1))
testdpqfun("norm2", c("d", "p", "q"), c(mean=1, sd=1))
testdpqfun("norm3", c("d", "p", "q"), c(mean=1, sd=1))

x <- rnorm(100)
fitdist(x, "norm") #ok
fitdist(x, "norm2", start=list(mean=1, sd=1)) #pnorm2 not defined
fitdist(x, "norm3", start=list(mean=1, sd=1)) #The dnorm3 function should return raise an error when names are incorrectly named

