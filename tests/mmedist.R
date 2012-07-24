library(fitdistrplus)



# (1) basic fit of a normal distribution with moment matching estimation
#

x1<-c(6.4, 13.3, 4.1, 1.3, 14.1, 10.6, 9.9, 9.6, 15.3, 22.1, 13.4, 
13.2, 8.4, 6.3, 8.9, 5.2, 10.9, 14.4)
mmedist(x1, "norm")

# (2) fit a discrete distribution (Poisson)
#

x2<-c(rep(4, 1), rep(2, 3), rep(1, 7), rep(0, 12))
mmedist(x2, "pois")



# (3) fit a finite-support distribution (beta)
#

x3<-c(0.80, 0.72, 0.88, 0.84, 0.38, 0.64, 0.69, 0.48, 0.73, 0.58, 0.81, 
0.83, 0.71, 0.75, 0.59)
mmedist(x3, "beta")


# (4) fit a Pareto distribution
#

if(any(installed.packages()[, "Package"] == "actuar"))
{
    require(actuar)
#simulate a sample
    x4 <- rpareto(1000, 6, 2)
	
#empirical raw moment
    memp <- function(x, order)
	ifelse(order == 1, mean(x), sum(x^order)/length(x))
	
	
#fit
    mmedist(x4, "pareto", order=c(1, 2), memp="memp", start=c(10, 10), 
			lower=1, upper=Inf)
	
}
