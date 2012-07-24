library(fitdistrplus)


# (1) for a fit of a normal distribution 
#

x1 <- c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
print(f1 <- fitdist(x1,"norm"))
gofstat(f1)
gofstat(f1,print.test=TRUE)

# (2) fit a discrete distribution (Poisson)
#

x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
print(f2<-fitdist(x2,"pois"))
g2 <- gofstat(f2,chisqbreaks=c(0,1),print.test=TRUE)
g2$chisqtable


# (3) comparison of fits of various distributions
#

x3<-rweibull(n=100,shape=2,scale=1)
gofstat(f3a<-fitdist(x3,"weibull"))
gofstat(f3b<-fitdist(x3,"gamma"))
gofstat(f3c<-fitdist(x3,"exp"))

# (4) Use of Chi-squared results in addition to
#     recommended statistics for continuous distributions
#

x4<-rweibull(n=100,shape=2,scale=1)
f4<-fitdist(x4,"weibull")
g4 <-gofstat(f4,meancount=10)
print(g4)

# (5) estimation of the standard deviation of a normal distribution 
# by maximum likelihood with the mean fixed at 10 using the argument fix.arg
#
f1b <- fitdist(x1,"norm",start=list(sd=5),fix.arg=list(mean=10))
gofstat(f1b)

