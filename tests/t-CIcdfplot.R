require("fitdistrplus")

nbboot <- 201
nbboot <- 10
ggplotEx <- requireNamespace("ggplot2", quietly = TRUE)

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
if(ggplotEx) CIcdfplot(b1, CI.level=95/100, CI.output = "probability", CI.fill="grey80", CI.col="black", plotstyle = "ggplot")
if(ggplotEx) CIcdfplot(b1, CI.level=95/100, CI.output = "quantile", datacol="blue", plotstyle = "ggplot")
  
par(mfrow=c(1,2))
CIcdfplot(b1, CI.level=85/100, CI.output = "probability")
CIcdfplot(b1, CI.level=85/100, CI.output = "quantile")
if(ggplotEx) CIcdfplot(b1, CI.level=85/100, CI.output = "probability", plotstyle = "ggplot")
if(ggplotEx) CIcdfplot(b1, CI.level=85/100, CI.output = "quantile", plotstyle = "ggplot")

par(mfrow=c(1,2))
CIcdfplot(b1, CI.level=90/100, CI.output = "probability")
CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="black", CI.type = "less", 
          CI.fill="grey85", verticals=TRUE, datacol="blue", do.points=FALSE)
if(ggplotEx) CIcdfplot(b1, CI.level=90/100, CI.output = "probability", plotstyle = "ggplot")
if(ggplotEx) CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="black", CI.type = "less", 
          CI.fill="grey85", verticals=TRUE, datacol="blue", do.points=FALSE, plotstyle = "ggplot")

par(mfrow=c(1,2))
CIcdfplot(b1, CI.level=90/100, CI.output = "probability", CI.type = "greater")
CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="black", CI.type = "greater", 
          CI.fill="grey90", datacol="blue", datapch=21)
if(ggplotEx) CIcdfplot(b1, CI.level=90/100, CI.output = "probability", CI.type = "greater", plotstyle = "ggplot")
if(ggplotEx) CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="black", CI.type = "greater", 
          CI.fill="grey90", datacol="blue", datapch=21, plotstyle = "ggplot")

par(mfrow=c(1,1))
CIcdfplot(b1, CI.level=90/100, CI.output = "probability", CI.col="black", CI.type = "less", CI.fill="grey90")
CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="black", CI.type = "less", CI.fill="grey90", 
          verticals=TRUE, datacol="blue", do.points=FALSE)
CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="grey90", CI.type = "less", CI.fill="grey90", 
          verticals=TRUE, datacol="blue", do.points=FALSE, CI.only=TRUE)
CIcdfplot(b1, CI.level=90/100, CI.output = "probability", CI.col="grey85", CI.type = "less", CI.fill="grey90", 
          CI.only = TRUE)
CIcdfplot(b1, CI.output = "probability", fitlty=3, fitlwd=4)

if(ggplotEx) CIcdfplot(b1, CI.level=90/100, CI.output = "probability", CI.col="black", CI.type = "less", CI.fill="grey90", plotstyle = "ggplot")
if(ggplotEx) CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="black", CI.type = "less", CI.fill="grey90", 
          verticals=TRUE, datacol="blue", do.points=FALSE, plotstyle = "ggplot")
if(ggplotEx) CIcdfplot(b1, CI.level=90/100, CI.output = "quantile", CI.col="grey90", CI.type = "less", CI.fill="grey90", 
          verticals=TRUE, datacol="blue", do.points=FALSE, CI.only=TRUE, plotstyle = "ggplot")
if(ggplotEx) CIcdfplot(b1, CI.level=90/100, CI.output = "probability", CI.col="grey85", CI.type = "less", CI.fill="grey90", 
          CI.only = TRUE, plotstyle = "ggplot")
if(ggplotEx) CIcdfplot(b1, CI.output = "probability", fitlty=3, fitlwd=4, plotstyle = "ggplot")

# (2) an example from ecotoxicology
# with censored data
#
data(salinity)
log10LC50 <-log10(salinity)
fln <- fitdistcens(log10LC50,"norm")
bln <- bootdistcens(fln, niter=nbboot)
(HC5ln <- quantile(bln,probs = 0.05))
CIcdfplot(bln, CI.output = "quantile", CI.fill = "lightblue", CI.col = "blue",
          xlab = "log10(LC50)",xlim=c(0.5,2),lines01 = TRUE)
if(ggplotEx) CIcdfplot(bln, CI.output = "quantile", CI.fill = "lightblue", CI.col = "blue",
          xlab = "log10(LC50)",xlim=c(0.5,2),lines01 = TRUE, plotstyle = "ggplot")

# zoom around the HC5 with CI on quantiles 
CIcdfplot(bln, CI.output = "quantile", CI.fill = "lightblue", CI.col = "blue",
          xlab = "log10(LC50)", lines01 = TRUE, xlim = c(0.8, 1.5), ylim = c(0, 0.1))
abline(h = 0.05, lty = 1)
if(ggplotEx) CIcdfplot(bln, CI.output = "quantile", CI.fill = "lightblue", CI.col = "blue",
          xlab = "log10(LC50)", lines01 = TRUE, xlim = c(0.8, 1.5), ylim = c(0, 0.1), plotstyle = "ggplot") + 
  ggplot2::geom_hline(yintercept = 0.05)

# zoom around the HC5 with CI on probabilities 
CIcdfplot(bln, CI.output = "probability", CI.fill = "lightblue", CI.col = "blue",
          xlab = "log10(LC50)", lines01 = TRUE, xlim = c(0.8, 1.5), ylim = c(0, 0.1))
abline(h = 0.05, lty = 1)
if(ggplotEx) CIcdfplot(bln, CI.output = "probability", CI.fill = "lightblue", CI.col = "blue",
          xlab = "log10(LC50)", lines01 = TRUE, xlim = c(0.8, 1.5), ylim = c(0, 0.1), plotstyle = "ggplot") + 
  ggplot2::geom_hline(yintercept = 0.05)

  
# (3) An example where the difference between "probability"
#     and "quantile" is clear on the plot
#

set.seed(123)
s3 <- rgamma(5, 3, 10)
f3 <- fitdist(s3, "norm")
b3 <- bootdist(f3, niter=nbboot, silent=TRUE)

par(mfrow=c(1,2))
CIcdfplot(b3, CI.level=90/100, CI.output = "probability")
CIcdfplot(b3, CI.level=90/100, CI.output = "quantile")

if(ggplotEx) CIcdfplot(b3, CI.level=90/100, CI.output = "probability", plotstyle = "ggplot")
if(ggplotEx) CIcdfplot(b3, CI.level=90/100, CI.output = "quantile", plotstyle = "ggplot")

#some ideas from http://edild.github.io/ssd/
