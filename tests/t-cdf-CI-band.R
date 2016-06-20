library(fitdistrplus)

nbboot <- 201
nbboot <- 10

# (1) Fit of a gamma distribution 
#

set.seed(123)
s1 <- rgamma(50, 3, 2)
f1 <- fitdist(s1, "gamma")
b1 <- bootdist(f1, niter=nbboot, silent=TRUE)

plot(b1)
quantile(b1)

par(mfrow=c(1,2))
CIcdfplot(b1, CI.level=90/100)
CIcdfplot(b1, CI.level=90/100, CI.col="grey80", CI.fill=TRUE, datacol="blue")


par(mfrow=c(1,2))
CIcdfplot(b1, CI.level=90/100, CI.type = "less")
CIcdfplot(b1, CI.level=90/100, CI.col="grey85", CI.type = "less", CI.fill=TRUE, 
        verticals=TRUE, datacol="blue", do.points=FALSE)

par(mfrow=c(1,2))
CIcdfplot(b1, CI.level=90/100, CI.type = "greater")
CIcdfplot(b1, CI.level=90/100, CI.col="grey90", CI.type = "greater", CI.fill=TRUE, 
        datacol="blue", datapch=21)

#some ideas from http://edild.github.io/ssd/
