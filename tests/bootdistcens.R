library(fitdistrplus)
nbboot <- 101
nbboot <- 11

# (1) Fit of a normal distribution to fluazinam data in log10
# followed by nonparametric bootstrap
#
data(fluazinam)
(d1 <-log10(fluazinam))
f1 <- fitdistcens(d1, "norm")
b1 <- bootdistcens(f1, niter = nbboot, silent=TRUE)
b1 <- bootdistcens(f1, niter = nbboot, silent=FALSE)


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

# (4) Comparison of fitdist and fitdistcens and bootdist and bootdistcens 
# for non censored data
x1<-c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
    13.2,8.4,6.3,8.9,5.2,10.9,14.4)
fx1<-fitdist(x1,"norm",method="mle")
cx1<-bootdist(fx1,bootmethod="nonparam", niter=nbboot)
xx1<-data.frame(left=x1,right=x1)
fxx1<-fitdistcens(xx1,"norm")
summary(fx1)
summary(fxx1)
cdfcomp(fx1)
cdfcompcens(fxx1)
cxx1<-bootdistcens(fxx1, niter=nbboot)
summary(cx1)
summary(cxx1)



# (5) fixing parameters
#
set.seed(1234)
x <- rexp(500, 5)
x <- data.frame(left=x, right=x+.1)

f1 <- fitdistcens(x, "gamma", fix.arg=list(shape=1.5))
b1 <- bootdistcens(f1, niter=nbboot)
plot(b1)

f1 <- fitdistcens(x, "gamma", fix.arg=function(x) list(shape=1.5))
b1 <- bootdistcens(f1, niter=nbboot)
plot(b1)

