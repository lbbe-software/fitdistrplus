library(fitdistrplus)
nbboot <- 100
nbboot <- 10
nsample <- 10

# (1) Description of a sample from a normal distribution
# with and without uncertainty on skewness and kurtosis estimated by bootstrap 
#
set.seed(1234)
x1 <- rnorm(nsample)
descdist(x1)
descdist(x1,boot=nbboot)

# (2) Description of a sample from a beta distribution
# with uncertainty on skewness and kurtosis estimated by bootstrap
# with changing of default colors 
#
descdist(rbeta(nsample,shape1=0.05,shape2=1),boot=nbboot,
obs.col="blue",boot.col="orange")

# (3) Description of a sample from a gamma distribution
# with uncertainty on skewness and kurtosis estimated by bootstrap
# without plotting 
#
descdist(rgamma(nsample,shape=2,rate=1),boot=nbboot,graph=FALSE)

# (3) Description of a sample from a Poisson distribution
# with uncertainty on skewness and kurtosis estimated by bootstrap 
#
descdist(rpois(nsample,lambda=2),discrete=TRUE,boot=nbboot)

# (4) Description of serving size data
# with uncertainty on skewness and kurtosis estimated by bootstrap 
#
data(groundbeef)
serving <- groundbeef$serving
descdist(serving, boot=nbboot)
