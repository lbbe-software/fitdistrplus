library(fitdistrplus)

# (1) Fit of a normal distribution followed by nonparametric bootstrap
#
d1<-data.frame(
left=c(1.73,1.51,0.77,1.96,1.96,-1.4,-1.4,NA,-0.11,0.55,
0.41,2.56,NA,-0.53,0.63,-1.4,-1.4,-1.4,NA,0.13),
right=c(1.73,1.51,0.77,1.96,1.96,0,-0.7,-1.4,-0.11,0.55,
0.41,2.56,-1.4,-0.53,0.63,0,-0.7,NA,-1.4,0.13))
f1<-fitdistcens(d1, "norm")
b1<-bootdistcens(f1)
b1
summary(b1)
plot(b1)

# (2) Fit of a gamma distribution followed by nonparametric bootstrap
#
d3<-data.frame(left=10^(d1$left),right=10^(d1$right))
f3 <- fitdistcens(d3,"gamma")
b3 <- bootdistcens(f3,niter=101)
summary(b3)
plot(b3)

# (3) Fit of a gamma distribution followed by nonparametric bootstrap
# with control of the optimization method
#
f3BFGS <- fitdistcens(d3,"gamma",optim.method="L-BFGS-B",lower=c(0,0))
b3BFGS <- bootdistcens(f3BFGS,niter=101)
summary(b3BFGS)
plot(b3BFGS)

# (4) Estimation of the standard deviation of a normal distribution 
# by maximum likelihood with the mean fixed at 0.1 using the argument fix.arg
# followed by nonparametric bootstrap
#
f1b <- fitdistcens(d1, "norm", start=list(sd=1.5),fix.arg=list(mean=0.1))
b1b<-bootdistcens(f1b,niter=101)
summary(b1b)
plot(b1b)
