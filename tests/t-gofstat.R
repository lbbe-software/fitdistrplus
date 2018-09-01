library(fitdistrplus)


# (1) fit of two distributions by maximum likelihood estimation
# to the serving size data
# and comparison of goodness-of-fit statistics
#

data(groundbeef)
serving <- groundbeef$serving
(fitg <- fitdist(serving, "gamma"))
gg <- gofstat(fitg)
(fitln <- fitdist(serving, "lnorm"))
gn <- gofstat(fitln)

gofstat(list(fitg, fitln))

# (2) fit of two discrete distributions to toxocara data
# and comparison of goodness-of-fit statistics
#

data(toxocara)
number <- toxocara$number

fitp <- fitdist(number, "pois")
summary(fitp)
plot(fitp)
gp <- gofstat(fitp)
gp

fitnb <- fitdist(number, "nbinom")
summary(fitnb)
plot(fitnb)
gnb <- gofstat(fitnb)
gnb

gofstat(list(fitp, fitnb))

attributes(gofstat(list(fitp, fitnb)))


# (3) Use of Chi-squared results in addition to
#     recommended statistics for continuous distributions
#

set.seed(1234)
x4 <- rweibull(n=10,shape=2,scale=1)
# fit of the good distribution
f4 <- fitdist(x4, "weibull")
g4  <- gofstat(f4, meancount=10)
print(g4)

# fit of a bad distribution
f4b <- fitdist(x4, "cauchy")
g4b  <- gofstat(f4b, meancount=10)
print(g4b)


# (4) estimation of the standard deviation of a normal distribution 
# by maximum likelihood with the mean fixed at 10 using the argument fix.arg
#
f1b <- fitdist(serving, "norm", start=list(sd=5), fix.arg=list(mean=10), lower=0)
gofstat(f1b)

