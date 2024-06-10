library(fitdistrplus)
nsample <- 10

# (1) basic fit of a normal distribution with moment matching estimation
#

set.seed(1234)
x1 <- rnorm(n=nsample)
mmedist(x1,"norm")

try(mmedist(x1,"norm", fix.arg=list(mean=0)))

# (2) fit a discrete distribution (Poisson)
#

set.seed(1234)
x2 <- rpois(n=nsample, lambda = 2)
mmedist(x2,"pois")

# (3) fit a finite-support distribution (beta)
#

set.seed(1234)
x3 <- rbeta(n=nsample, shape1=5, shape2=10)
mmedist(x3,"beta")



# (4) fit a Pareto distribution
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
  mmedist(x4, "pareto", order=c(1, 2), memp=memp, start=c(shape=10, scale=10), 
          lower=1, upper=Inf)
  mmedist(x4, "pareto", order=1, memp=memp, start=list(shape=10), fix.arg=list(scale=1.5), 
          lower=2, upper=Inf)
  mmedist(x4, "pareto", order=1, memp=memp, start=function(x) list(shape=10), fix.arg=list(scale=1.5), 
          lower=2, upper=Inf)
  mmedist(x4, "pareto", order=1, memp=memp, start=list(shape=10), fix.arg=function(x) list(scale=1.5), 
          lower=2, upper=Inf)
  
  #weights
  memp2 <- function(x, order, weights) sum(x^order * weights)/sum(weights)
  w <- rep(1, length(x4))
  w[x4 < 1] <- 2
  mmedist(x4, "pareto", order=c(1, 2), memp=memp2, weights=w,
          start=list(shape=10, scale=10), lower=1, upper=Inf)
  
  #fit
  data(danishuni)
  fparedanishMME <- mmedist(danishuni$Loss, "pareto", order=1:2, 
                            memp=memp, start=c(shape=10, scale=10), lower=2+1e-6, upper=Inf)
  c(theo = mpareto(1, fparedanishMME$estimate[1], fparedanishMME$estimate[2]),
    emp = memp(danishuni$Loss, 1))	
  c(theo = mpareto(2, fparedanishMME$estimate[1], fparedanishMME$estimate[2]),
    emp = memp(danishuni$Loss, 2))
  
}


# (5) fit a lognormal distribution
#
x3 <- rlnorm(n=nsample)
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
  jensen=as.numeric(exp(2*f1$estimate["meanlog"]+f1$estimate["sdlog"]^2)*(exp(f1$estimate["sdlog"]^2)-1)), 
  emp=var(x3)*(n-1)/n)


# (6) test error messages
#

mnorm3 <- dnorm3 <- function(x, a)
  "NA"
x <- rnorm(nsample)

#should get a one-line error 
res <- mmedist(x, "norm3", start=list(a=1), order=1, memp=function(x, order) mean(x))
#as in 
attr(try(log("a"), silent=TRUE), "condition")


try(mmedist(c(serving, "a"), "gamma"))
try(mmedist(c(serving, NA), "gamma"))
try(mmedist(c(serving, Inf), "gamma"))
try(mmedist(c(serving, -Inf), "gamma"))
try(mmedist(c(serving, NaN), "gamma"))


# (7) fit of a normal distribution with weighted moment matching estimation
#

n <- nsample
w <- c(rep(1, n/2), rep(10, n/2))
mmedist(x1, "norm", weights=w)$estimate

#check
sum(w*x1)/sum(w)
fitdistrplus:::wtdmean(x1, w)
sum(w*(x1-sum(w*x1)/sum(w))^2)/sum(w)
fitdistrplus:::wtdvar(x1, w)


mmedist(exp(x1), "lnorm", weights=w)$estimate

#test non integer weights
try(mmedist(x1, "norm", weights=rep(1/3, length(x1))))
try(mmedist(1:10, "pois", weights=c(rep(1, 9), 1.001), start=list(lambda=mean(x))))
try(mmedist(1:10, "pois", weights=c(rep(1, 9), 1.0000001), start=list(lambda=mean(x))))


# (8) fit of a neg binom distribution with weighted moment matching estimation
#

x4 <- rnbinom(nsample, 5, 1/2)
table(x4)
w <- rep(1, length(x4))
w[x4 > 5] <- 2

mmedist(x4, "nbinom", weights=w)$estimate
mmedist(x4, "nbinom", weights=NULL)$estimate

# (9) relevant example for zero modified geometric distribution
#

memp1 <- function(x, order) mean(x^order)
memp2 <- function(x, order, weights) sum(x^order * weights)/sum(weights)


rzmgeom <- function(n, p1, p2)
{
  u <- rbinom(n, 1, 1-p1) #prob to get zero is p1
  u[u != 0] <- rgeom(sum(u != 0), p2)+1
  u
}
dzmgeom <- function(x, p1, p2)
{
  p1 * (x == 0) + (1-p1)*dgeom(x-1, p2)
}
mcumulgeom <- function(order, prob) #E[N(N-1)..(N-o+1)]
{
  (1-prob)^order * factorial(order) / prob^order
}
mgeom <- function(order, prob)
{
  if(order == 0)
    return(1)
  m1 <- mcumulgeom(1, prob)
  if(order == 1)
    return(m1)
  m2 <- mcumulgeom(2, prob) + m1
  if(order == 2)
    return(m2)
  m3 <- mcumulgeom(3, prob) + 3*m2 + 2*m1
  if(order == 3)
    return(m3)
  m4 <- mcumulgeom(4, prob) + 6*m3 + 11*m2 - 6*m1
  if(order == 4)
    return(m4)
  else
    stop("not yet implemented")
}
if(FALSE)
{
  for(j in 1:4)
    print(c(memp1(rgeom(1e4, 1/6), j), mgeom(j, 1/6)))
}

mzmgeom <- function(order, p1, p2) #raw moment
{
  momj <- sapply(0:order, function(j) mgeom(j, p2))
  cj <- choose(order, 0:order)
  (1-p1) * sum( cj * momj )
}
if(FALSE)
{
  for(j in 1:4)
    print(c(memp1(rzmgeom(1e4, 1/3, 1/6), j), mzmgeom(j, 1/3, 1/6)))
  
}


x5 <- rzmgeom(nsample, 1/3, 1/6)
w <- rep(1, length(x5))
w[x5 > 5] <- 2

mmedist(x5, "zmgeom", order=1:2, memp=memp1, start=list(p1=mean(x5 == 0), p2=1/mean(x5[x5 > 0])), 
        lower=0.01, upper=0.99)$estimate
mmedist(x5, "zmgeom", order=1:2, memp=memp2, start=list(p1=mean(x5 == 0), p2=1/mean(x5[x5 > 0])),  
        weights=w)$estimate


mmedist(x5, "zmgeom", order=1:2, memp=memp1, start=list(p1=mean(x5 == 0), p2=1/mean(x5[x5 > 0])), 
        lower=0.01, upper=0.99)$loglik
mmedist(x5, "zmgeom", order=1:2, memp=memp2, start=list(p1=mean(x5 == 0), p2=1/mean(x5[x5 > 0])),  
        weights=w)$loglik

# (10) bounds
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
  mmedist(x4, "pareto", order=c(1, 2), memp=memp, start=c(shape=10, scale=10), lower=1, upper=Inf, optim.method = "L-BFGS-B") #L-BFGS-B via optim
  mmedist(x4, "pareto", order=c(1, 2), memp=memp, start=c(shape=10, scale=10), lower=1, upper=Inf, optim.method = "Nelder") #Nelder Mead via constrOptim
}
