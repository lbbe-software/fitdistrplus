# startargdefault function returns initial values of parameters generally using moments or quantiles

# INPUTS 
#x : data vector or matrix
#distr : the distribution name

# OUTPUTS
# a named list or raises an error 
startargdefault <- function(x, distr)
{
  if(distr %in% c("fpareto", "pareto4", "pareto3", "pareto2", "pareto1", "pareto",
                  "genpareto", "trbeta", "burr", "llogis", "paralogis", "invburr",
                  "invpareto",  "invparalogis"))
  {
    start <- startarg_fellerpareto_family(x, distr)
  }else if(distr %in% c("trgamma", "gamma", "weibull", "exp"))
  {
    start <- startarg_transgamma_family(x, distr)
  }else if(distr %in% c("invtrgamma", "invgamma", "invweibull", "invexp"))
  {
    start <- startarg_invtransgamma_family(x, distr)
  }else if (distr == "norm") {
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
  }else if (distr == "nbinom") {
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    size <- ifelse(v > m, m^2/(v - m), 100)
    start <- list(size = size, mu = m) 
  }else if (distr == "geom" ) {
    m <- mean(x)
    prob <- ifelse(m>0, 1/(1+m), 1)
    start <- list(prob=prob)        
  }else if (distr == "beta") {
    if (any(x < 0) | any(x > 1)) 
      stop("values must be in [0-1] to fit a beta distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    aux <- m*(1-m)/v - 1
    start <- list(shape1=m*aux, shape2=(1-m)*aux)
  }else if (distr == "logis") {
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    start <- list(location=m, scale=sqrt(3*v)/pi)
  }else if (distr == "cauchy") {
    start <- list(location=median(x), scale=IQR(x)/2)
  }else if (distr == "unif"){
    start <- list(min=0, max=1)
  }else if (distr == "lgamma")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-gamma distribution")
    #p228 of Klugmann and Hogg (1984)
    m1 <- mean(log(x))
    m2 <- mean(log(x)^2)
    alpha <- m1^2/(m2-m1^2)
    lambda <- m1/(m2-m1^2)
    start <- list(shapelog=alpha, ratelog=lambda)
  }else
    stop(paste0("Unknown starting values for distribution ", distr, "."))
  
  return(start)
} 

startarg_fellerpareto_family <- function(x, distr)
{
  distlist <- c("fpareto", "pareto4", "pareto3", "pareto2", "pareto1", "pareto",
                "genpareto", "trbeta", "burr", "llogis", "paralogis", "invburr",
                "invpareto",  "invparalogis")
  stopifnot(length(distr) == 1)
  stopifnot(distr %in% distlist)
  
  gammarel <- function(alpha, q1, q3)
  {
    T1 <- ((4/3)^(1/alpha) - 1) / (4^(1/alpha) - 1)
    T2 <- (q1 - muhat)/(q3 - muhat)
    log(T1)/log(T2)
  }
  thetarel <- function(alpha, gamma, q1, q3)
  {
    (q1 - q3)/( ((4/3)^(1/alpha) - 1)^(1/gamma) - (4^(1/alpha) - 1)^(1/gamma) )
  }  
  alpharel <- function(x, mu, theta, gamma)
  {
    y <- ((x-mu)/theta)^(gamma)
    1/mean(log(1+y))
  }
  eps <- 5/100
  
  if(distr == "pareto4")
  {
    muhat <- min(x, na.rm=TRUE)
    if(muhat < 0)
      muhat <- muhat*(1+eps)
    else
      muhat <- muhat*(1-eps)
    alphastar <- 2
    q1 <- quantile(x, probs=1/4, na.rm=TRUE)
    q3 <- quantile(x, probs=3/4, na.rm=TRUE)
    gammahat <- as.numeric(gammarel(alphastar, q1, q3))
    thetahat <- as.numeric(thetarel(alphastar, gammahat, q1, q3))
    alphahat <- alpharel(x, muhat, thetahat, gammahat)
    
    start <- list(min=muhat, shape1=alphahat, shape2=gammahat, scale=thetahat)
  }else if(distr == "pareto3")
  {
    muhat <- min(x, na.rm=TRUE)
    if(muhat < 0)
      muhat <- muhat*(1+eps)
    else
      muhat <- muhat*(1-eps)
    alphastar <- 1 #true value in that case
    q1 <- quantile(x, probs=1/4, na.rm=TRUE)
    q3 <- quantile(x, probs=3/4, na.rm=TRUE)
    gammahat <- as.numeric(gammarel(alphastar, q1, q3))
    thetahat <- as.numeric(thetarel(alphastar, gammahat, q1, q3))
    
    start <- list(min=muhat, shape=gammahat, scale=thetahat)
    
  }else if(distr == "pareto2")
  {
    muhat <- min(x, na.rm=TRUE)
    if(muhat < 0)
      muhat <- muhat*(1+eps)
    else
      muhat <- muhat*(1-eps)
    gammastar <- 1 #true value
    alphastar <- 2
    q1 <- quantile(x, probs=1/4, na.rm=TRUE)
    q3 <- quantile(x, probs=3/4, na.rm=TRUE)
    thetahat <- as.numeric(thetarel(alphastar, gammastar, q1, q3))
    alphahat <- alpharel(x, muhat, thetahat, gammastar)
    
    start <- list(min=muhat, shape=alphahat, scale=thetahat)
  }else if (distr == "pareto1")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a Pareto distribution")
    muhat <- min(x, na.rm=TRUE)
    if(muhat < 0)
      muhat <- muhat*(1+eps)
    else
      muhat <- muhat*(1-eps)
    alphat <- 1/mean(log(x/muhat), na.rm = TRUE)
    
    start <- list(shape=alphat, min=muhat)
  }else if(distr == "pareto")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a Pareto distribution")
    m1 <- mean(x, na.rm = TRUE)
    m2 <- mean(x^2, na.rm = TRUE)
    alphastar <- max(2*(m2-m1^2)/(m2-2*m1^2), 2)
    mustar <- 0 #true value
    gammastar <- 1 #true value
    q1 <- quantile(x, probs=1/4, na.rm=TRUE)
    q3 <- quantile(x, probs=3/4, na.rm=TRUE)
    thetahat <- as.numeric(thetarel(alphastar, gammastar, q1, q3))
    alphahat <- alpharel(x, mustar, thetahat, gammastar)
    
    start <- list(shape=alphahat, scale=thetahat)
  }else if (distr == "llogis")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-logistic  distribution")
    q25 <- as.numeric(quantile(x, 0.25, na.rm=TRUE))
    q75 <- as.numeric(quantile(x, 0.75, na.rm=TRUE))
    shape <- 2*log(3)/(log(q75)-log(q25))
    scale <- exp(log(q75)+log(q25))/2
    start <- list(shape=shape, scale=scale)
  }else
    stop("wrong distr")
  return(start)
}

startarg_transgamma_family <- function(x, distr)
{
  distlist <- c("trgamma", "gamma", "weibull", "exp")
  stopifnot(length(distr) == 1)
  stopifnot(distr %in% distlist)
  
  if (distr == "trgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit a trans-gamma distribution")
    #same as gamma with shape2=tau=1
    m <- mean(x, na.rm = TRUE)
    v <- var(x, na.rm = TRUE)
    alphahat <- m^2/v
    thetahat <- v/m
    tauhat <- digamma(alphahat) / mean(log(x/thetahat))
    
    start <- list(shape1=alphahat, shape2=tauhat, scale=thetahat)
  }else if (distr == "gamma") 
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a gamma distribution")
    #same as gamma with shape2=tau=1
    m <- mean(x, na.rm = TRUE)
    v <- var(x, na.rm = TRUE)
    alphahat <- m^2/v
    thetahat <- v/m
    
    start <- list(shape=alphahat, scale=thetahat)
  }else if (distr == "weibull") 
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a Weibull distribution")
    q25 <- as.numeric(quantile(x, 0.25, na.rm = TRUE))
    q75 <- as.numeric(quantile(x, 0.75, na.rm = TRUE))
    if(q25 < q75) #check to avoid division by zero
    {
      q50 <- median(x, na.rm = TRUE)
      tauhat <- (log(-log(1-1/4)) - log(-log(1-3/4))) / (log(q25) - log(q75))
      thetahat <- q50/(-log(1-1/2))^(tauhat)
    }else
    {
      m <- mean(log(x), na.rm = TRUE)
      v <- var(log(x), na.rm = TRUE)
      tauhat <- 1.2/sqrt(v)
      thetahat <- exp(m + 0.572/tauhat)
    }
    
    start <- list(shape = tauhat, scale = thetahat)
  }else if (distr == "exp") 
  {
    if (any(x < 0)) 
      stop("values must be positive to fit an exponential distribution")
    start <- list(rate=1/mean(x))
  }else
    stop("wrong distr")
  return(start)
}

startarg_invtransgamma_family <- function(x, distr)
{
  distlist <- c("invtrgamma", "invgamma", "invweibull", "invexp")
  stopifnot(length(distr) == 1)
  stopifnot(distr %in% distlist)
  
  
 if (distr == "invtrgamma") {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse trans-gamma distribution")
    start <- startarg_transgamma_family(1/x, "trgamma")
  }else if (distr == "invgamma")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse gamma  distribution")
    #http://en.wikipedia.org/wiki/Inverse-gamma_distribution
    m1 <- mean(x, na.rm = TRUE)
    m2 <- mean(x^2, na.rm = TRUE)
    shape <- (2*m2-m1^2)/(m2-m1^2)
    scale <- m1*m2/(m2-m1^2)
    start <- list(shape=shape, scale=scale)
  }else if (distr == "invweibull")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse Weibull distribution")
    g <- log(log(4))/(log(log(4/3)))
    p25 <- as.numeric(quantile(x, 0.25))
    p75 <- as.numeric(quantile(x, 0.75))
    shape <- exp((g*log(p75)-log(p25))/(g-1))
    scale <-log(log(4))/(log(shape)-log(p25))
    start <- list(shape=shape, scale=max(scale, 1e-9))
  }else if(distr == "invexp")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse exponential distribution")
    start <- startarg_transgamma_family(1/x, "exp")
  }else
    stop("wrong distr")
  return(start)
}