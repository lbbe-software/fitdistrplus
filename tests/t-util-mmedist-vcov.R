require(fitdistrplus)


#### analytical example for dfapprox ####
if(FALSE)
{
  f <- function(x)
    sin(x[1]) + cos(x[2])
  fprime <- function(x, i)
  {
    if(i == 1)
      return(cos(x[1]))
    -sin(x[2])
  }
  
  x <- c(pi/4, pi/4)
  fitdistrplus:::dfapprox(f, 1, x) - fprime(x, 1)
  
  
  dfapprox_1 <- function(x)
    sapply(x, function(y) fitdistrplus:::dfapprox(f, 1, y))
  
  f <- function(x)
    sin(x)
  fprime <- function(x)
    cos(x)
  
  x <- seq(0, pi, pi/10)
  dfapprox_1(x) - fprime(x)
  
  curve(dfapprox_1(x), from=-pi/2, to=pi/2)
  curve(fprime(x), add=TRUE, col="green")
  
}

#### Gamma example ####
if(FALSE)
{
  memp  <-  function(x, order) mean(x^order)
  
  #true value, see Ibragimov & Has'minskii (1981)
  J <- function(alpha, beta)
  {
    cbind(c(1/beta, (1+2*alpha)/beta^2), 
          c(-alpha/beta^2, -2*alpha*(alpha+1)/beta^3))
  }
  
  J(3, 1/2)
  fitdistrplus:::mme.Jhat(c(3, 1/2), NULL, order=1:2, mgamma)
  
  
  diff <- function(x)
    max(abs(J(x[1], x[2]) - fitdistrplus:::mme.Jhat(x, NULL, order=1:2, mgamma)))
  
  x <- seq(1/4, 3/2, by=1/4)
  
  #shape param
  sapply(x, function(x) diff(c(x, 1/2)))
  #rate param
  sapply(x, function(x) diff(c(1/2, x)))
  
  
  #true value, see Ibragimov & Has'minskii (1981)
  A <- function(alpha, beta)
  {
    cbind(c(alpha/beta^2, 2*alpha*(alpha+1)/beta^3),
          c(2*alpha*(alpha+1)/beta^3, alpha*(4*alpha+6)*(alpha+1)/beta^4))
  }
  
  n <- 1e6
  x <- rgamma(n, 3, 1/2)
  A(3, 1/2)
  fitdistrplus:::mme.Ahat(c(3, 1/2), fix.arg=NULL, order=1:2, obs=x, mdistnam=mgamma, memp, weights=NULL)
  
  
  vcovgamma <- function(alpha, beta)
  {
    Jinv <- solve(J(alpha, beta))
    Jinv %*% A(alpha, beta) %*% t(Jinv)
  }
  
  #two param, no fix
  fitdistrplus:::mme.vcov(c(3, 1/2), fix.arg=NULL, order=1:2, obs=x, mdistnam=mgamma, memp, weights=NULL)
  vcovgamma(3, 1/2)
  
  #one param, one fix
  fitdistrplus:::mme.vcov(3, fix.arg=list(rate=1/2), order=1, obs=x, mdistnam=mgamma, memp, weights=NULL)

  #cancel out a rown and a line in the inverse matrix
  Jinv <- 1/J(3, 1/2)
  Jinv[2,1] <- Jinv[, 2] <- 0
  Jinv %*% A(3, 1/2) %*% t(Jinv)

}


