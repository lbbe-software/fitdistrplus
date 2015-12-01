library(fitdistrplus)


#(1) beta distribution
#

n <- 100
x <- rbeta(n, 3, 3/4)



par(mfrow=c(1, 3))
llsurface(plot.min=c(0.1, 0.1), plot.max=c(7, 3), plot.arg=c("shape1", "shape2"),
                         plot.np=50, obs=x, distr="beta", plot.type="persp", theta=60, phi=30)
llsurface(plot.min=c(0.1, 0.1), plot.max=c(7, 3), plot.arg=c("shape1", "shape2"),
                         plot.np=50, obs=x, distr="beta", plot.type="contour")
points(3, 3/4, pch="+", col="red")
llsurface(plot.min=c(0.1, 0.1), plot.max=c(7, 3), plot.arg=c("shape1", "shape2"),
                         plot.np=50, obs=x, distr="beta", plot.type="image")



par(mfrow=c(1,2))
llcurve(plot.min=0.1, plot.max=7, plot.arg="shape1", fix.arg=list(shape2=3/4), plot.np=50, 
                      obs=x, distr="beta")
llcurve(plot.min=0.1, plot.max=2, plot.arg="shape2", fix.arg=list(shape1=3), plot.np=50, 
                      obs=x, distr="beta")
