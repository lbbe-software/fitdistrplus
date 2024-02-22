library(actuar)
library(fitdistrplus)

n <- 1e2
U <- rgamma(n, 1/2)
V <- rgamma(n, 3)
X <- rfpareto(n, 0, 1/2, 3/2, 5/2, scale=2)

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

thetarel2 <- function(x, mu)
{
  y <- 2*x*(x-mu)/(1+(x-mu)^2)
  z <- x/(x-mu)
  sqrt(mean(y)/(mean(z)-1))
}


thetahat <- thetarel2(X, min(X)*.95)

Y <- (X-min(X)*.95)/thetahat
Z <- Y/(1+Y)
Y <- (X-3)/10
Z <- Y/(1+Y)

cdfcomp(fitdist(Z^2, "beta"))
cdfcomp(fitdist(Z, "beta"))
coef(fitdist(Z, "beta"))
coef(fitdist(Z^2, "beta"))

tauhat <- coef(fitdist(Z, "beta"))[1]
alphahat <- coef(fitdist(Z, "beta"))[2]

gammarel <- function(x, mu, theta, alpha, tau)
{
  y <- log(1+(x-mu)/theta)
  z <- log((x-mu)/theta)
  1+(digamma(tau) - digamma(alpha+tau) + mean(y))/mean(z)
}

gammarel(X, min(X)*.9, thetahat, alphahat, tauhat)
