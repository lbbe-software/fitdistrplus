library(actuar)
library(fitdistrplus)

n <- 1e2
U <- rgamma(n, 1/2)
V <- rgamma(n, 3)
X <- rfpareto(n, 0, 1/2, 3/2, 5/2, 1)

hist(log(U/V))
hist(-log(1-U/V))

Z <- log(U/V)
cdfcomp(fitdist(Z, "norm"))


hist(log(X))
hist(X)
hist(X/(1+X))
hist(log((1+X)/X))

Y <- X-min(X)*.95
Z <- Y/(1+Y)
cdfcomp(fitdist(Z, "beta"))
coef(fitdist(Z, "beta"))


X <- rfpareto(n, 3, 1/2, 3/2, 5/2, 10)

thetarel <- function(x, alpha, gamma)
{
  q1 <- quantile(x, probs=1/4)
  q3 <- quantile(x, probs=3/4)
  
  ( ((3/4)^(-1/alpha) - 1)^(1/gamma) - ((1/4)^(-1/alpha) - 1)^(1/gamma) )/(q1 - q3)
}  

thetahat <- thetarel(X, 2, 2)

Y <- (X-min(X)*.95)/thetahat
Z <- Y/(1+Y)
Y <- (X-min(X)*.95)/sd(X)
Z <- Y/(1+Y)
Y <- (X-3)/10
Z <- Y/(1+Y)

cdfcomp(fitdist(Z, "beta"))
coef(fitdist(Z, "beta"))

