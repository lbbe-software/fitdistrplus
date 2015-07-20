#x : data vector or matrix
#distr : the distribution name
start.arg.default <- function(x, distr)
{
  if (distr == "norm") {
    n <- length(x)
    sd0 <- sqrt((n - 1)/n) * sd(x)
    mx <- mean(x)
    start <- list(mean=mx, sd=sd0)
  }else if (distr == "lnorm") {
    if (any(x <= 0)) 
      stop("values must be positive to fit a lognormal distribution")
    n <- length(x)
    lx <- log(x)
    sd0 <- sqrt((n - 1)/n) * sd(lx)
    ml <- mean(lx)
    start <- list(meanlog=ml, sdlog=sd0)
  }else if (distr == "pois") {
    start <- list(lambda=mean(x))
  }else if (distr == "exp") {
    start <- list(rate=1/mean(x))
  }else if (distr == "gamma") {
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    start <- list(shape=m^2/v, rate=m/v)
  }else if (distr == "nbinom") {
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    size <- if (v > m) m^2/(v - m)
    else 100
    start <- list(size = size, mu = m) 
  }else if (distr == "geom" ) {
    m <- mean(x)
    prob <- if (m>0) 1/(1+m)
    else 1
    start <- list(prob=prob)        
  }else if (distr == "beta") {
    if (any(x < 0) | any(x > 1)) 
      stop("values must be in [0-1] to fit a beta distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    aux <- m*(1-m)/v - 1
    start <- list(shape1=m*aux, shape2=(1-m)*aux)
  }else if (distr == "weibull") {
    m <- mean(log(x))
    v <- var(log(x))
    shape <- 1.2/sqrt(v)
    scale <- exp(m + 0.572/shape)
    start <- list(shape = shape, scale = scale)
  }else if (distr == "logis") {
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    start <- list(location=m, scale=sqrt(3*v)/pi)
  }else if (distr == "cauchy") {
    start <- list(location=median(x), scale=IQR(x)/2)
  }else if (distr == "unif"){
    start <- list(min=0, max=1)
  }else if (distr == "invgamma")
  {
    #http://en.wikipedia.org/wiki/Inverse-gamma_distribution
    m1 <- mean(x)
    m2 <- mean(x^2)
    shape <- (2*m2-m1^2)/(m2-m1^2)
    scale <- m1*m2/(m2-m1^2)
    start <- list(shape=shape, scale=scale)
  }else if (distr == "llogis")
  {
    p25 <- quantile(x, 0.25)
    p75 <- quantile(x, 0.75)
    shape <- 2*log(3)/(log(p75)-log(p25))
    scale <- exp(log(p75)+log(p25))/2
    start <- list(shape=shape, scale=scale)
  }else if (distr == "invweibull")
  {
    g <- log(log(4))/(log(log(4/3)))
    p25 <- quantile(x, 0.25)
    p75 <- quantile(x, 0.75)
    shape <- exp((g*log(p75)-log(p25))/(g-1))
    scale <-log(log(4))/(log(shape)-log(p75))
    start <- list(shape=shape, scale=scale)
  }else if (distr == "pareto1")
  {
    #http://www.math.umt.edu/gideon/pareto.pdf
    x1 <- min(x)
    m1 <- mean(x)
    n <- length(x)
    shape <- (n*m1-x1)/(n*(m1-x1))
    min <- x1*(n*shape - 1)/(n*shape)
    start <- list(shape=shape, min=min)
  }else if (distr == "pareto2")
  {
    m1 <- mean(x)
    m2 <- mean(x^2)
    scale <- (m1*m2)/(m2-2*m1^2)
    shape <- 2*(m2-m1^2)/(m2-2*m1^2)
    start <- list(shape=shape, scale=scale)
  }else
    stop(paste0("Unknown starting values for distribution ", distr, "."))
  
  return(start)
} 