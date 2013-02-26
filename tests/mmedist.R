library(fitdistrplus)

# (1) basic fit of a normal distribution with moment matching estimation
#

set.seed(1234)
x1 <- rnorm(n=100)
mmedist(x1,"norm")

# (2) fit a discrete distribution (Poisson)
#

set.seed(1234)
x2 <- rpois(n=30,lambda = 2)
mmedist(x2,"pois")

# (3) fit a finite-support distribution (beta)
#

set.seed(1234)
x3 <- rbeta(n=100,shape1=5, shape2=10)
mmedist(x3,"beta")



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
    mmedist(x4, "pareto", order=c(1, 2), memp="memp", start=c(shape=10, scale=10), 
			lower=1, upper=Inf)
			
#fit
data(danishuni)
fparedanishMME <- fitdist(danishuni$Loss, "pareto", method="mme", order=1:2, 
      memp="memp", start=c(shape=10, scale=10), lower=2+1e-6, upper=Inf, 
      control=list(trace=1))
c(theo = mpareto(1, fparedanishMME$estimate[1], fparedanishMME$estimate[2]),
emp = memp(danishuni$Loss, 1))	
c(theo = mpareto(2, fparedanishMME$estimate[1], fparedanishMME$estimate[2]),
emp = memp(danishuni$Loss, 2))
	
}


# (5) fit a lognormal distribution
#

f1 <- mledist(x3, "lnorm") #previously mmedist was the same as mledist
f2 <- mmedist(x3, "lnorm")
n <- length(x3)
s2 <- log(1+var(x3)/mean(x3)^2*(n-1)/n)
mu <- log(mean(x3)) - s2/2
cbind(c(mu, s2), f2$estimate)


c(truestim=exp(mu+s2/2), 
	jensen=as.numeric(exp(f1$estimate["meanlog"]+f1$estimate["sdlog"]^2/2)), 
	emp=mean(x3))

c(truestim=exp(2*mu+s2)*(exp(s2)-1), 
	jensen=as.numeric(exp(f1$estimate["meanlog"]+f1$estimate["sdlog"]^2/2)*(exp(f1$estimate["sdlog"]^2)-1)), 
	emp=var(x3)*(n-1)/n)


