require("fitdistrplus")

# (1) Plot of an empirical censored distribution (censored data) as a CDF
# using the default Turnbull method
#
data(smokedfish)
plotdistcens(smokedfish)
plotdistcens(data.frame(right=smokedfish$right, left=smokedfish$left))
d1 <- as.data.frame(log10(smokedfish))
plotdistcens(d1)

#test on first arg
try(plotdistcens(list(left=smokedfish$left, right=smokedfish$right)))
try(plotdistcens(cbind(left=smokedfish$left, right=smokedfish$right)))
d2 <- data.frame(left=smokedfish$right, right=smokedfish$left)
try(plotdistcens(d2))

# (2) Add the CDF of a normal distribution and QQ and PP plots
#
plotdistcens(smokedfish,"lnorm", para=list(meanlog=-3.6,sdlog=3.5))
plotdistcens(d1,"norm", para=list(mean=-1.6,sd=1.5))

# (3) Various plots of the same empirical distribution 
#
# default Wang plot
plotdistcens(d1, NPMLE = TRUE, NPMLE.method = "Wang")
plotdistcens(d1, NPMLE = TRUE, NPMLE.method = "Wang", lwd = 3, main = "Wang ECDF plot")
# Turnbull plot
plotdistcens(d1, NPMLE = TRUE, NPMLE.method = "Turnbull", col = "red", 
             main = "Turnbull ECDF plot")
plotdistcens(d1,Turnbull = TRUE) # deprecated way to do it
# Turnbull plot with confidence intervals
plotdistcens(d1,NPMLE = TRUE, NPMLE.method = "Turnbull", Turnbull.confint = TRUE)
plotdistcens(d1,Turnbull = TRUE,Turnbull.confint = TRUE) # deprecated way to do it
# with intervals and points
plotdistcens(d1,NPMLE = FALSE)
plotdistcens(d1,NPMLE = FALSE, col = "red", lwd = 2)
plotdistcens(d1,rightNA=3, NPMLE = FALSE)
plotdistcens(d1,rightNA=3, Turnbull = FALSE) # deprecated way to do it
# with intervals and points
# defining a minimum value for left censored values
plotdistcens(d1,leftNA=-3, NPMLE = FALSE)

# (4) Goodness-of-fit plots for the same dataset after logarithmic transformation
#   with a lognormal distribution, successively using the three proposed methods
#
d3 <- smokedfish
plotdistcens(d3,"lnorm",para=list(meanlog=-3.6,sdlog=3.5), main = "Wang plot")
plotdistcens(d3,"lnorm",para=list(meanlog=-3.6,sdlog=3.5), 
             NPMLE.method = "Turnbull", main = "Turnbull plot")
plotdistcens(d3,"lnorm",para=list(meanlog=-3.6,sdlog=3.5),
             NPMLE = FALSE, leftNA=0, main = "Plot of ordered intervals")

# Test with the salinity data set
#
data(salinity)
log10LC50 <-log10(salinity)
plotdistcens(log10LC50)
plotdistcens(log10LC50, NPMLE.method = "Turnbull")
plotdistcens(log10LC50, NPMLE = FALSE)
fn <- fitdistcens(log10LC50,"norm")
fl <- fitdistcens(log10LC50,"logis")
plot(fn)
plot(fl)
