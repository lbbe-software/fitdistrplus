library(actuar)
library(fitdistrplus)

n <- 1e2
x <- rpareto4(n, theopar["mu"], theopar["alpha"], theopar["gamma"], 1/theopar["theta"])

theopar <- c("mu"=1, "alpha"=.75, "gamma"=3, "theta"=4)

#check constant
gammarel1 <- function(p, min, shape1, shape2, rate)
{
  p1 <- p
  p2 <- 1-p
  q1 <- qpareto4(p1,  min, shape1, shape2, rate)
  q2 <- qpareto4(p2,  min, shape1, shape2, rate)
  T1 <- ((1-p1)^(-1/shape1) - 1) / ((1-p2)^(-1/shape1) - 1)
  T2 <- (q1 - min)/(q2 - min)
  #print(T1)
  #print(T2)
  log(T1)/log(T2)
}
gammarel1(1/4, theopar["mu"], theopar["alpha"], theopar["gamma"], 1/theopar["theta"])

curve(gammarel1(x, theopar["mu"], theopar["alpha"], theopar["gamma"], 1/theopar["theta"]),
      from=.1, to= .49)
abline(h=theopar["gamma"], lty=2)

gammarel3 <- function(p, min, shape1, shape2, rate)
{
  p1 <- p
  p2 <- 1-p
  q1 <- quantile(x, p1)
  q2 <- quantile(x, p2)
  T1 <- ((1-p1)^(-1/shape1) - 1) / ((1-p2)^(-1/shape1) - 1)
  T2 <- (q1 - min)/(q2 - min)
  #print(T1)
  #print(T2)
  log(T1)/log(T2)
}

curve(gammarel3(x, theopar["mu"], theopar["alpha"], theopar["gamma"], 1/theopar["theta"]),
      from=.1, to= .49)
abline(h=theopar["gamma"], lty=2)


thetarel1 <- function(p, min, shape1, shape2, rate)
{
  p1 <- p
  p2 <- 1-p
  q1 <- qpareto4(p1,  min, shape1, shape2, rate)
  q2 <- qpareto4(p2,  min, shape1, shape2, rate)
  T1 <- ((1-p1)^(-1/shape1) - 1)^(1/shape2) - ((1-p2)^(-1/shape1) - 1)^(1/shape2)
  T1/(q1-q2)
}

curve(thetarel1(x, theopar["mu"], theopar["alpha"], theopar["gamma"], 1/theopar["theta"]),
      from=.1, to= .49)
abline(h=theopar["theta"], lty=2)


cbind("theo"=mpareto4(1:5, theopar["mu"], theopar["alpha"], theopar["gamma"], 1/theopar["theta"]),
      "emp"=sapply(1:5, function(j) mean(x^j)))

coef(fitdist(x/(1+x), "beta"))
coef(fitdist((x-1)/x, "beta"))

q1 <- quantile(x, probs=1/4)
q3 <- quantile(x, probs=3/4)

muhat <- min(x)*0.95

gammarel2 <- function(x)
{
  T1 <- ((3/4)^(-1/x) - 1) / ((1/4)^(-1/x) - 1)
  T2 <- (q1 - muhat)/(q3 - muhat)
  log(T1)/log(T2)
}
thetarel <- function(alpha, gamma)
{
  (q1 - q3)/( ((3/4)^(-1/alpha) - 1)^(1/gamma) - ((1/4)^(-1/alpha) - 1)^(1/gamma) )
}  
alpharel <- function(x, mu, theta, gamma)
{
  y <- ((x-mu)/theta)^(gamma)
  1/mean(log(1+y))
}


curve(gammarel2(x), from=0.5, to=5, ylim=c(0, 6))
abline(h=1/2)
abline(h=theopar["gamma"])

gammahat1 <- gammarel2(1)
thetahat1 <- thetarel(1, gammahat1)
alphahat1 <- alpharel(x, muhat, thetahat1, gammahat1)
gammahat1bis <- gammarel2(alphahat1)
thetahat1bis <- thetarel(alphahat1, gammahat1bis)

gammahat2 <- gammarel2(2)
thetahat2 <- thetarel(2, gammahat2)
alphahat2 <- alpharel(x, muhat, thetahat2, gammahat2)

gammahat3 <- gammarel2(3)
thetahat3 <- thetarel(3, gammahat3)
alphahat3 <- alpharel(x, muhat, thetahat3, gammahat3)

cbind(theopar, 
      "alpha=1"=c(muhat, alphahat1, gammahat1, thetahat1),
      "alpha=1, bis"=c(muhat, alphahat1, gammahat1bis, thetahat1bis),
      "alpha=2"=c(muhat, alphahat2, gammahat2, thetahat2),
      "alpha=3"=c(muhat, alphahat3, gammahat3, thetahat3)
      )


theopar <- c("mu"=1, "alpha"=2, "gamma"=3, "theta"=4)


geomean <- function(n)
{
  x <- rpareto4(n, theopar["mu"], theopar["alpha"], theopar["gamma"], 1/theopar["theta"])
  log(prod(x^(1/n)) / min(x))
}


gn <- replicate(1e4, geomean(n))
hist(gn)
mean(gn)

x <- rpareto1(n, theopar["alpha"], theopar["mu"])
gn <- replicate(n, prod(rpareto1(n, theopar["alpha"], theopar["mu"])^(1/n)))
hist(gn)
abline(v=exp(1/theopar["alpha"]))
mean(gn)


