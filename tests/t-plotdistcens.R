library(fitdistrplus)

# (1) Plot of an empirical censored distribution (censored data) as a CDF
# using the default Turnbull method
#
data(smokedfish)
d1 <- as.data.frame(log10(smokedfish))
plotdistcens(d1)

# (2) Add the CDF of a normal distribution and QQ and PP plots
#
plotdistcens(d1,"norm",para=list(mean=-1.6,sd=1.5))

# (3) Various plots of the same empirical distribution 
#
# default Wang plot
plotdistcens(d1, NPMLE = TRUE, NPMLE.method = "Wang")
# Turnbull plot
plotdistcens(d1, NPMLE = TRUE, NPMLE.method = "Turnbull")
plotdistcens(d1,Turnbull = TRUE) # deprecated way to do it
# Turnbull plot with confidence intervals
plotdistcens(d1,NPMLE = TRUE, NPMLE.method = "Turnbull", Turnbull.confint = TRUE)
plotdistcens(d1,Turnbull = TRUE,Turnbull.confint = TRUE) # deprecated way to do it
# with intervals and points
plotdistcens(d1,rightNA=3, NPMLE = FALSE)
plotdistcens(d1,rightNA=3, Turnbull = FALSE) # deprecated way to do it
# with intervals and points
# defining a minimum value for left censored values
plotdistcens(d1,leftNA=-3, NPMLE = FALSE)

# (4) Plot of the CDF of the same dataset after logarithmic transformation
#   with a lognormal distribution, successively using the two proposed methods
#
d3<-data.frame(left=10^(d1$left),right=10^(d1$right))
plotdistcens(d3,"lnorm",para=list(meanlog=0.27,sdlog=3.3))
plotdistcens(d3,"lnorm",para=list(meanlog=0.27,sdlog=3.3),NPMLE = FALSE, leftNA=0)

