library(fitdistrplus)

# (1) non-censored data (continuous)
#

s1 <- NULL
s2 <- list("mean"=2)
s0 <- list("mean"=2, "sd"=1)
s3 <- function(x)
  list("mean"=1.01*mean(x)) 
s4 <- function(x)
  list("mean"=1.01*mean(x), "sd"=sd(x)) 

f1 <- NULL  
f2 <- list("sd"=3)
f3 <- function(x) 
  list("sd"=1.01*sd(x))  

x <- rnorm(1000)

#redefine normal distribution for better check
dnorm2 <- dnorm
pnorm2 <- pnorm
qnorm2 <- qnorm
rnorm2 <- rnorm
mnorm2 <- function(order, mean, sd)
  ifelse(order == 1, mean, sd^2)
memp  <-  function(x, order) mean(x^order)

# both NULL
mf1 <- mledist(x, "norm", start=s1, fix.arg=f1)
mf1 <- mmedist(x, "norm", start=s1, fix.arg=f1)
mf1 <- qmedist(x, "norm", start=s1, fix.arg=f1, probs=1:2/3)
mf1 <- mgedist(x, "norm", start=s1, fix.arg=f1)

fit1 <- fitdist(x, "norm", start=s1, fix.arg=f1)
boot1 <- bootdist(fit1, niter=10)

# both named list
mf1 <- mledist(x, "norm2", start=s2, fix.arg=f2)
mf1 <- mmedist(x, "norm2", start=s2, fix.arg=f2, order=1, memp=memp)
mf1 <- qmedist(x, "norm2", start=s2, fix.arg=f2, probs=1/3)
mf1 <- mgedist(x, "norm2", start=s2, fix.arg=f2)

fit1 <- fitdist(x, "norm2", start=s2, fix.arg=f2)
boot1 <- bootdist(fit1, niter=10)

# named list and NULL
mf1 <- mledist(x, "norm2", start=s0, fix.arg=f1)
mf1 <- mmedist(x, "norm2", start=s0, fix.arg=f1, order=1:2, memp=memp)
mf1 <- qmedist(x, "norm2", start=s0, fix.arg=f1, probs=1:2/3)
mf1 <- mgedist(x, "norm2", start=s0, fix.arg=f1)

fit1 <- fitdist(x, "norm2", start=s0, fix.arg=f1)
boot1 <- bootdist(fit1, niter=10)

# NULL and named list 
mf1 <- mledist(x, "norm", start=s1, fix.arg=f2)
mf1 <- qmedist(x, "norm", start=s1, fix.arg=f2, probs=1/3)
mf1 <- mgedist(x, "norm", start=s1, fix.arg=f2)

fit1 <- fitdist(x, "norm", start=s1, fix.arg=f2)
boot1 <- bootdist(fit1, niter=10)

# both function
mf1 <- mledist(x, "norm2", start=s3, fix.arg=f3)
mf1 <- mmedist(x, "norm2", start=s3, fix.arg=f3, order=1, memp=memp)
mf1 <- qmedist(x, "norm2", start=s3, fix.arg=f3, probs=1/3)
mf1 <- mgedist(x, "norm2", start=s3, fix.arg=f3)

fit1 <- fitdist(x, "norm2", start=s3, fix.arg=f3)
boot1 <- bootdist(fit1, niter=10)

# function and NULL
mf1 <- mledist(x, "norm2", start=s4, fix.arg=f1)
mf1 <- mmedist(x, "norm2", start=s4, fix.arg=f1, order=1:2, memp=memp)
mf1 <- qmedist(x, "norm2", start=s4, fix.arg=f1, probs=1:2/3)
mf1 <- mgedist(x, "norm2", start=s4, fix.arg=f1)

fit1 <- fitdist(x, "norm2", start=s4, fix.arg=f1)
boot1 <- bootdist(fit1, niter=10)

# NULL and function
mf1 <- mledist(x, "norm", start=s1, fix.arg=f3)
mf1 <- qmedist(x, "norm", start=s1, fix.arg=f3, probs=1/3)
mf1 <- mgedist(x, "norm", start=s1, fix.arg=f3)

fit1 <- fitdist(x, "norm", start=s1, fix.arg=f3)
boot1 <- bootdist(fit1, niter=10)

# should raise error for too less parameters
try(mgedist(x, "norm", start=s2, fix.arg=f1))
try(fitdist(x, "norm", start=s2, fix.arg=f1))
# should raise error for too much parameters
try(mgedist(x, "norm", start=s0, fix.arg=f2))
try(fitdist(x, "norm", start=s0, fix.arg=f2))
# should raise error for NA value
try(mgedist(x, "norm", start=s1, fix.arg=list(sd=NA)))
try(fitdist(x, "norm", start=list(sd=NA)))
# should raise error for inconsistent parameter
try(mgedist(x, "norm", start=function(x) list("toto"=1)))
try(fitdist(x, "norm", fix=list(toto=2)))

# (2) censored data
#

data(salinity)
log10LC50 <-log10(salinity)

s1 <- NULL
s2 <- list("mean"=2)
s0 <- list("mean"=2, "sd"=1)
s3 <- function(x) list("mean"=mean(x))

f1 <- NULL  
f2 <- list("sd"=3)
f3 <- function(x) list("sd"=sd(x))

fitdistcens(log10LC50, "norm", start=s1, fix.arg = f1)
fitdistcens(log10LC50, "norm", start=s1, fix.arg = f2)
fitdistcens(log10LC50, "norm", start=s2, fix.arg = f2)
fitdistcens(log10LC50, "norm", start=s0, fix.arg = f1)
fitdistcens(log10LC50, "norm", start=s3, fix.arg = f2)
fitdistcens(log10LC50, "norm", start=s3, fix.arg = f3)
fitdistcens(log10LC50, "norm", start=s1, fix.arg = f3)


fit1 <- fitdistcens(log10LC50, "norm", start=s1, fix.arg = f1)
boot1 <- bootdistcens(fit1, niter = 10)

fit1 <- fitdistcens(log10LC50, "norm", start=s3, fix.arg = f2)
boot1 <- bootdistcens(fit1, niter = 10)

fit1 <- fitdistcens(log10LC50, "norm", start=s2, fix.arg = f3)
boot1 <- bootdistcens(fit1, niter = 10)

# (3) non-censored data (discrete)
#

n <- 200
trueval <- c("size"=10, "prob"=3/4, "mu"=10/3)
x <- rnbinom(n, trueval["size"], trueval["prob"])

mledist(x, "nbinom")
fitdist(x, "nbinom")


# (4) non-censored data (continuous) external distributions
#

data("endosulfan")
ATV <-endosulfan$ATV
fendo.ln <- fitdist(ATV, "lnorm")
fendo.g <- fitdist(ATV, "gamma", start=list(shape=2, scale=1), lower=0)
require("actuar")
fendo.ll <- fitdist(ATV, "llogis", start = list(shape = 1, scale = 500))
fendo.P <- fitdist(ATV, "pareto", start = list(shape = 1, scale = 500))
fendo.B <- fitdist(ATV, "burr", start = list(shape1 = 0.3, shape2 = 1, 
                                             rate = 1))
