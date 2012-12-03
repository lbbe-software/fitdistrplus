library(fitdistrplus)
nbboot <- 101

# (1) Fit of a normal distribution to fluazinam data in log10
# followed by nonparametric bootstrap
#
data(fluazinam)
(d1 <-log10(fluazinam))
f1 <- fitdistcens(d1, "norm")
b1 <- bootdistcens(f1, niter = nbboot)
b1
summary(b1)
plot(b1)

# (2) Estimation of the mean of the normal distribution 
# by maximum likelihood with the standard deviation fixed at 1 
# using the argument fix.arg
# followed by nonparametric bootstrap with less iterations
#
f1b <- fitdistcens(d1, "norm", start = list(mean = 1), fix.arg = list(sd = 1))
b1b <- bootdistcens(f1b, niter = nbboot)
summary(b1b)
plot(b1b)

# (3) Estimation of the standard deviation of a normal distribution 
# by maximum likelihood with the mean fixed at 0.1 using the argument fix.arg
# followed by nonparametric bootstrap
#
f1b <- fitdistcens(d1, "norm", start=list(sd=1.5), fix.arg=list(mean=0.1))
b1b <- bootdistcens(f1b, niter=nbboot)
summary(b1b)
plot(b1b)
