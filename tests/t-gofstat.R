require("fitdistrplus")

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

# (5) Use on a small data set (less than 10 observations)
# no pb identified
#

set.seed(1234)
x5a <- rweibull(n=4,shape=2,scale=1)
f5a <- fitdist(x5a, "weibull")
(g5a  <- gofstat(f5a))

x5b <- rpois(n = 4, lambda = 1)
f5b <- fitdist(x5b, "pois")
(g5b <- gofstat(f5b))

nsample <- 500
nsample <- 10
visualize <- FALSE # TRUE for manual tests with visualization of results
set.seed(1234)

# (6) censored dataset
#
data(fluazinam)
log10EC50 <-log10(fluazinam)

# call fitdistcens with a 'custom' optimization function
fit.1 <- fitdistcens(log10EC50, "logis")
fit.2 <- fitdistcens(log10EC50, "norm")

gofstat(fit.1, fitnames = "logistic")
gofstat(list(fit.1, fit.2), fitnames = c("logistic", "normal"))
