require(fitdistrplus)
nsample <- 1e6
nsample <- 10

#### (1) Gamma example ####

truetheta <- c("alpha"=3, "beta"=1/2)
x <- rgamma(nsample, truetheta["alpha"], truetheta["beta"])
f1 <- mmedist(x, "gamma", order=1:2, calcvcov = TRUE)
f1$vcov
if(FALSE)
{
  memp  <-  function(x, order) mean(x^order)
  require(actuar)
  fitdistrplus:::mme.vcov(as.numeric(truetheta), fix.arg=NULL, order=1:2, obs=x, mdistnam=mgamma, memp, weights=NULL)
}


# (2) fit a Pareto distribution
#

if(any(installed.packages()[, "Package"] == "actuar"))
{
  require(actuar)
  #simulate a sample
  x4 <- rpareto(nsample, 6, 2)
  
  #empirical raw moment
  memp <- function(x, order)
    mean(x^order)
  
  #fit
  res <- mmedist(x4, "pareto", order=c(1, 2), memp=memp, start=c(shape=10, scale=10), 
          lower=1, upper=Inf, calcvcov = TRUE)
}