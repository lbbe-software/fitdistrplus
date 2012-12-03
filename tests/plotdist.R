library(fitdistrplus)

# (1) Plot of an empirical distribution with changing 
# of default line types for CDF and colors
#
set.seed(1234)
x1 <- rnorm(n=30)
plotdist(x1)
plotdist(x1,col="red",type="b",pch=16)
plotdist(x1,type="s")

# (2) Plot of a discrete distribution against data
#
set.seed(1234)
x2 <- rpois(n=30,lambda = 2)
plotdist(x2,discrete=TRUE)
plotdist(x2,"pois",para=list(lambda=mean(x2)))
plotdist(x2,"pois",para=list(lambda=mean(x2)),col="red",lwd="2")

# (3) Plot of a continuous distribution against data
#
xn <- rnorm(n=100,mean=10,sd=5)
plotdist(xn,"norm",para=list(mean=mean(xn),sd=sd(xn)))
plotdist(xn,"norm",para=list(mean=mean(xn),sd=sd(xn)),pch=16,col="green")

# (4) Plot of serving size data
#
data(groundbeef)
plotdist(groundbeef$serving,type="s")

