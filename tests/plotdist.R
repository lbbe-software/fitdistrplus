library(fitdistrplus)



# (1) Plot of an empirical distribution with changing 
# of default line types for CDF and colors
#
x1<-c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
plotdist(x1)
plotdist(x1,col="red",type="b",pch=16)
plotdist(x1,type="s")

# (2) Plot of a discrete distribution against data
#
x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
plotdist(x2,discrete=TRUE)
plotdist(x2,"pois",para=list(lambda=mean(x2)))
plotdist(x2,"pois",para=list(lambda=mean(x2)),col="red",lwd="2")

# (3) Plot of a continuous distribution against data
#
xn<-rnorm(n=100,mean=10,sd=5)
plotdist(xn,"norm",para=list(mean=mean(xn),sd=sd(xn)))
plotdist(xn,"norm",para=list(mean=mean(xn),sd=sd(xn)),pch=16,col="green")

# (4) Plot of serving size data
#
data(groundbeef)
plotdist(groundbeef$serving)
