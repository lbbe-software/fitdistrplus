library(fitdistrplus)



# (1) Plot of an empirical censored distribution (censored data) as a CDF
# using the default Turnbull method
#
d1<-data.frame(
left=c(1.73,1.51,0.77,1.96,1.96,-1.4,-1.4,NA,-0.11,0.55,
0.41,2.56,NA,-0.53,0.63,-1.4,-1.4,-1.4,NA,0.13),
right=c(1.73,1.51,0.77,1.96,1.96,0,-0.7,-1.4,-0.11,0.55,
0.41,2.56,-1.4,-0.53,0.63,0,-0.7,NA,-1.4,0.13))
plotdistcens(d1)
plotdistcens(d1,col="red")

# (2) Add the CDF of a normal distribution 
#
plotdistcens(d1,"norm",para=list(mean=0.12,sd=1.4))

# (3) Basic plot of the same empirical distribution with intervals and points
# defining a realistic maximum value for right censored values
# in the second plot
#
plotdistcens(d1,Turnbull = FALSE)
plotdistcens(d1,rightNA=3, Turnbull = FALSE)

# (4) Plot of the CDF of the same dataset after logarithmic transformation
#   with a lognormal distribution, successively using the two proposed methods
#
d3<-data.frame(left=10^(d1$left),right=10^(d1$right))
plotdistcens(d3,"lnorm",para=list(meanlog=0.27,sdlog=3.3))
plotdistcens(d3,"lnorm",para=list(meanlog=0.27,sdlog=3.3),Turnbull = FALSE, leftNA=0)

