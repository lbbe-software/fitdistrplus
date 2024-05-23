require(fitdistrplus)
nsample <- 10

#### (1) Gamma example ####

require(actuar)

truetheta <- c("alpha"=3, "beta"=1/2)
x <- rgamma(nsample, truetheta["alpha"], truetheta["beta"])
f1 <- fitdist(x, "gamma", method="mme", order=1:2)
summary(f1)

# (4) fit a Pareto distribution
#

#if(any(installed.packages()[, "Package"] == "actuar"))
{
  require(actuar)
  #simulate a sample
  x4 <- rpareto(nsample, 6, 2)
  
  #empirical raw moment
  memp <- function(x, order)
    mean(x^order)
  
  #fit
  mmedist(x4, "pareto", order=c(1, 2), memp=memp, start=c(shape=10, scale=10), 
          lower=1, upper=Inf)
}