library(fitdistrplus)

nbboot <- 201
nbboot <- 10

# (1) Fit of a gamma distribution 
#

set.seed(123)
s1 <- rgamma(20, 3, 2)
f1 <- fitdist(s1, "gamma")
b1 <- bootdist(f1, niter=nbboot, silent=TRUE)

plot(b1)
quantile(b1)

par(mfrow=c(1,2))
CIcdfplot(b1, CI.level=95/100, CI.output = "probability", CI.fill="grey80", CI.col="black")
CIcdfplot(b1, CI.level=95/100, CI.output = "quantile", datacol="blue")


par(mfrow=c(1,2))
CIcdfplot(b1, CI.level=90/100, CI.output = "probability")
CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="black", CI.type = "less", 
          CI.fill="grey85", verticals=TRUE, datacol="blue", do.points=FALSE)

par(mfrow=c(1,2))
CIcdfplot(b1, CI.level=90/100, CI.output = "probability", CI.type = "greater")
CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="black", CI.type = "greater", 
          CI.fill="grey90", datacol="blue", datapch=21)


par(mfrow=c(1,1))
CIcdfplot(b1, CI.level=90/100, CI.output = "probability", CI.col="black", CI.type = "less", CI.fill="grey90")
CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="black", CI.type = "less", CI.fill="grey90", 
          verticals=TRUE, datacol="blue", do.points=FALSE)
CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="grey90", CI.type = "less", CI.fill="grey90", 
          verticals=TRUE, datacol="blue", do.points=FALSE, CI.only=TRUE)
CIcdfplot(b1, CI.level=90/100, CI.output = "probability", CI.col="grey85", CI.type = "less", CI.fill="grey90", 
          CI.only = TRUE)

# (2) an example from ecotoxicology
# with censored data
#
data(salinity)
log10LC50 <-log10(salinity)
fln <- fitdistcens(log10LC50,"norm")
bln <- bootdistcens(fln, niter=101)
(HC5ln <- quantile(bln,probs = 0.05))
CIcdfplot(bln, CI.output = "quantile", CI.fill = "lightblue", CI.col = "blue",
          xlab = "log10(LC50)",xlim=c(0.5,2),lines01 = TRUE)

# zoom around the HC5 with CI on quantiles 
# visual problem near 0
CIcdfplot(bln, CI.output = "quantile", CI.fill = "lightblue", CI.col = "blue",
          xlab = "log10(LC50)", lines01 = TRUE, xlim = c(0.8, 1.5), ylim = c(0, 0.1))
abline(h = 0.05, lty = 1)

# zoom around the HC5 with CI on probabilities 
CIcdfplot(bln, CI.output = "probability", CI.fill = "lightblue", CI.col = "blue",
          xlab = "log10(LC50)", lines01 = TRUE, xlim = c(0.8, 1.5), ylim = c(0, 0.1))
abline(h = 0.05, lty = 1)


#some ideas from http://edild.github.io/ssd/
