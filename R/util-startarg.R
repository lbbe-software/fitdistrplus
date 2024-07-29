#############################################################################
#   Copyright (c) 2024 Christophe Dutang, Marie Laure Delignette-Muller                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################
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
  }else if(distr %in% c("lgamma", "gumbel", "invgauss", "genbeta"))
  {
    start <- startarg_othercontinuous_actuar_family(x, distr)
  }else if(distr %in% c("ztpois", "ztbinom", "ztnbinom", "ztgeom",
                         "logarithmic", "zmpois", "zmbinom", "zmnbinom",
                         "zmgeom", "zmlogarithmic", "poisinvgauss")) 
  {
    start <- startarg_discrete_actuar_family(x, distr)
  }else if (distr == "norm") {
    n <- length(x)
    sd0 <- sqrt((n - 1)/n) * sd(x, na.rm = TRUE)
    mx <- mean(x, na.rm = TRUE)
    start <- list(mean=mx, sd=sd0)
  }else if (distr == "lnorm") {
    if (any(x <= 0)) 
      stop("values must be positive to fit a lognormal distribution")
    n <- length(x)
    lx <- log(x)
    sd0 <- sqrt((n - 1)/n) * sd(lx, na.rm = TRUE)
    ml <- mean(lx, na.rm = TRUE)
    start <- list(meanlog=ml, sdlog=sd0)
  }else if (distr == "pois") {
    start <- list(lambda=mean(x))
  }else if (distr == "binom") {
    m <- mean(x, na.rm = TRUE)
    v <- var(x, na.rm = TRUE)
    prob <- ifelse(v < m, 1-v/m, 1)
    size <- ceiling(max(x, m/prob, na.rm = TRUE))
    start <- list(size = size, prob = prob)
  }else if (distr == "nbinom") {
    n <- length(x)
    m <- mean(x, na.rm = TRUE)
    v <- (n - 1)/n*var(x, na.rm = TRUE)
    size <- ifelse(v > m, m^2/(v - m), 100)
    start <- list(size = size, mu = m) 
  }else if (distr == "geom" ) {
    m <- mean(x, na.rm = TRUE)
    prob <- ifelse(m>0, 1/(1+m), 1)
    start <- list(prob=prob)        
  }else if (distr == "beta") {
    if (any(x < 0) | any(x > 1)) 
      stop("values must be in [0-1] to fit a beta distribution")
    n <- length(x)
    m <- mean(x, na.rm = TRUE)
    v <- (n - 1)/n*var(x, na.rm = TRUE)
    aux <- m*(1-m)/v - 1
    start <- list(shape1=m*aux, shape2=(1-m)*aux)
  }else if (distr == "logis") {
    n <- length(x)
    m <- mean(x, na.rm = TRUE)
    v <- (n - 1)/n*var(x, na.rm = TRUE)
    start <- list(location=m, scale=sqrt(3*v)/pi)
  }else if (distr == "cauchy") {
    q2 <- median(x, na.rm = TRUE)
    start <- list(location=q2, scale=IQR(x, na.rm = TRUE)/2)
  }else if (distr == "unif"){
    start <- list(min=min(x), max=max(x))
  }else
    stop(paste0("Unknown starting values for distribution ", distr, "."))
  
  return(start)
} 

startarg_discrete_actuar_family <- function(x, distr)
{
  distlist <- c("ztpois", "ztbinom", "ztnbinom", "ztgeom",
                "logarithmic", "zmpois", "zmbinom", "zmnbinom",
                "zmgeom", "zmlogarithmic", "poisinvgauss")
  stopifnot(length(distr) == 1)
  stopifnot(distr %in% distlist)
  
  if (distr == "logarithmic")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a logarithmic distribution")
    m1 <- mean(x, na.rm = TRUE)
    start <- list(prob = 1-1/m1)
  }else if (distr == "ztpois")
  {
    if (any(x < 1)) 
      stop("values must be greater than 1 to fit a ZT-Poisson distribution")
    MMEpois <- startargdefault(x-1, "pois")
    start <- list(lambda = MMEpois$lambda)
  }else if (distr == "ztnbinom")
  {
    if (any(x < 1)) 
      stop("values must be greater than 1 to fit a ZT-negative binomial distribution")
    MMEnbinom <- startargdefault(x-1, "nbinom")
    prob <- MMEnbinom$size/(MMEnbinom$size + MMEnbinom$mu)
    start <- list(size = MMEnbinom$size, prob = prob) 
  }else if (distr == "ztgeom")
  {
    if (any(x < 1)) 
      stop("values must be greater than 1 to fit a ZT-geometric distribution")
    MMEgeom <- startargdefault(x-1, "geom")
    start <- list(prob = MMEgeom$prob) 
  }else if (distr == "ztbinom")
  {
    if (any(x < 1)) 
      stop("values must be greater than 1 to fit a ZT-binomial distribution")
    MMEbinom <- startargdefault(x-1, "binom")
    start <- list(size = MMEbinom$size, prob = MMEbinom$prob) 
  }else if (distr == "zmlogarithmic")
  {
    p0 <- mean(x == 0, na.rm = TRUE)
    MMElog <- startarg_discrete_actuar_family(x, "logarithmic")
    
    start <- list(p0 = p0, prob = MMElog$prob*(1-p0)) 
  }else if (distr == "zmpois")
  {
    p0 <- mean(x == 0, na.rm = TRUE)
    MMEpois <- startargdefault(x, "pois")
    
    start <- list(p0 = p0, lambda = MMEpois$lambda*(1-p0)) 
  }else if (distr == "zmnbinom")
  {
    p0 <- mean(x == 0, na.rm = TRUE)
    MMEnbinom <- startargdefault(x, "nbinom")
    prob <- MMEnbinom$size/(MMEnbinom$size + MMEnbinom$mu)
    start <- list(p0 = p0, size = MMEnbinom$size, prob = prob*(1-p0)) 
  }else if (distr == "zmgeom")
  {
    p0 <- mean(x == 0, na.rm = TRUE)
    MMEgeom <- startargdefault(x, "geom")
    start <- list(p0=p0, prob = MMEgeom$prob*(1-p0)) 
  }else if (distr == "zmbinom")
  {
    p0 <- mean(x == 0, na.rm = TRUE)
    MMEbinom <- startargdefault(x, "binom")
    start <- list(p0=p0, size = MMEbinom$size, prob = MMEbinom$prob*(1-p0)) 
  }else if(distr == "poisinvgauss") 
  {
    m <- mean(x == 0, na.rm = TRUE)
    v <- var(x == 0, na.rm = TRUE)
    phihat <- ifelse(v > m, (v-m)/m^3, v/m^3)
    start <- list(mean = m, dispersion = phihat)
    
  }else
    stop(paste0("Unknown starting values for distribution ", distr, "."))
  
  return(start)
}

startarg_othercontinuous_actuar_family <- function(x, distr)
{
  distlist <- c("lgamma", "gumbel", "invgauss", "genbeta")
  stopifnot(length(distr) == 1)
  stopifnot(distr %in% distlist)
  if (distr == "lgamma")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-gamma distribution")
    #p228 of Klugmann and Hogg (1984)
    m1 <- mean(log(x), na.rm = TRUE)
    m2 <- mean(log(x)^2, na.rm = TRUE)
    alpha <- m1^2/(m2-m1^2)
    lambda <- m1/(m2-m1^2)
    start <- list(shapelog=alpha, ratelog=lambda)
  }else if(distr == "gumbel")
  {
    q1 <- as.numeric(quantile(x, probs=1/4, na.rm=TRUE))
    q2 <- median(x, na.rm = TRUE)
    q3 <- as.numeric(quantile(x, probs=3/4, na.rm=TRUE))
    thetahat <- (q1 - q3)/(log(log(4/3)) - log(log(4)))
    alphahat <- thetahat*log(log(2)) + q2
    start <- list(scale = thetahat, alpha=alphahat)
  }else if(distr == "invgauss")
  {
    m1 <- mean(x, na.rm = TRUE)
    v <- var(x, na.rm = TRUE)
    start <- list(mean=m1, dispersion = v/m1^3)
  }else if(distr == "genbeta")
  {
    m1 <- mean(x, na.rm = TRUE)
    m2 <- mean(x^2, na.rm = TRUE)
    taustar <- 3
    thetahat <- m2/m1
    if(thetahat > 1)
    {
      thetahat <- thetahat*max(x, na.rm = TRUE)
      y <- (x/thetahat)^taustar
    }else
    {
      thetahat <- 1/thetahat*max(x, na.rm = TRUE)
      y <- (x/thetahat)^taustar
    }
    betaMME <- startargdefault(y, "beta")
    
    start <- list(shape1=betaMME$shape1, shape2 = betaMME$shape2,
                  shape3=taustar, scale=thetahat)
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
  #Pareto 4
  gammarelP4 <- function(alpha, q1, q3)
  {
    T1 <- ((4/3)^(1/alpha) - 1) / (4^(1/alpha) - 1)
    T2 <- (q1 - muhat)/(q3 - muhat)
    log(T1)/log(T2)
  }
  thetarelP4 <- function(alpha, gamma, q1, q3)
  {
    (q1 - q3)/( ((4/3)^(1/alpha) - 1)^(1/gamma) - (4^(1/alpha) - 1)^(1/gamma) )
  }  
  alpharelP4 <- function(x, mu, theta, gamma)
  {
    y <- ((x-mu)/theta)^(gamma)
    1/mean(log(1+y), na.rm = TRUE)
  }
  #Feller Pareto
  thetarelFP <- function(x, mu)
  {
    y <- 2*x*(x-mu)/(1+(x-mu)^2)
    z <- x/(x-mu)
    sqrt(mean(y, na.rm = TRUE)/(mean(z, na.rm = TRUE)-1))
  }
  gammarelFP <- function(x, mu, theta, alpha, tau)
  {
    y <- log(1+(x-mu)/theta)
    z <- log((x-mu)/theta)
    1+(digamma(tau) - digamma(alpha+tau) + mean(y, na.rm = TRUE))/mean(z, na.rm = TRUE)
  }
  #transformed beta
  thetarelTRB <- function(x)
  {
    y <- 2*x^2/(1+x^2)
    sqrt(mean(y, na.rm = TRUE))
  }
  #generalized Pareto
  thetarelGP <- function(x)
  {
    mean(x/(1+x), na.rm = TRUE)
  }
  #Burr
  thetarelBurr <- function(x)
  {
    y <- 2*x^2/(1+x^2)
    sqrt(2/3/mean(y, na.rm = TRUE))
  }
  #paralogistic
  alpharelPL <- function(x)
  {
    y <- log(x^2/(1+x^2))
    z <- log(x)*x^2/(1+x^2)
    (mean(y, na.rm = TRUE) - mean(z, na.rm = TRUE))/(mean(z, na.rm = TRUE) - 2)
  }
  thetarelPL <- function(x, alpha)
  {
    y <- x^(alpha+1)/(1+x^alpha)
    (alpha+1)*mean(y, na.rm = TRUE)
  }
  
  eps <- 5/100
  murel <- function(x)
  {
    muhat <- min(x, na.rm=TRUE)
    if(muhat < 0)
      muhat <- muhat*(1+eps)
    else
      muhat <- muhat*(1-eps)
    muhat
  }  
  
  if(distr == "fpareto")
  {
    muhat <- murel(x)
    thetahat <- thetarelFP(x, muhat)
    
    y <- (x-muhat)/thetahat
    z <- y/(1+y)
    betaMME <- startargdefault(z, "beta")
    tauhat <- betaMME$shape1
    alphahat <- betaMME$shape2
    gammahat <- gammarelFP(x, muhat, thetahat, alphahat, tauhat)
      
    start <- list(min=muhat, shape1=alphahat, shape2=gammahat, shape3=tauhat, scale=thetahat)
  }else if(distr == "trbeta")
  {
    thetahat <- thetarelTRB(x)
    
    y <- (x)/thetahat
    z <- y/(1+y)
    betaMME <- startargdefault(z, "beta")
    tauhat <- betaMME$shape1
    alphahat <- betaMME$shape2
    mustar <- 0 #true value
    gammahat <- gammarelFP(x, mustar, thetahat, alphahat, tauhat)
    
    start <- list(shape1=alphahat, shape2=gammahat, shape3=tauhat, scale=thetahat)
  }else if(distr == "genpareto")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a Burr distribution")
    thetahat <- thetarelGP(x)
    y <- (x)/thetahat
    z <- y/(1+y)
    betaMME <- startargdefault(z, "beta")
    tauhat <- betaMME$shape1
    alphahat <- betaMME$shape2
    
    start <- list(shape1=alphahat, shape2=tauhat, scale=thetahat)
  }else if(distr == "burr")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a Burr distribution")
    thetahat <- thetarelBurr(x)
    taustar <- 1 #true value
    alphastar <- 2
    mustar <- 0 #true value
    gammahat <- gammarelFP(x, mustar, thetahat, alphastar, taustar)
    q2 <- median(x, na.rm=TRUE)
    alphahat <- log(2)/log(1+(q2/thetahat)^gammahat)
    
    start <- list(shape1=alphahat, shape2=gammahat, scale=thetahat)
  }else if (distr == "invburr")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a inverse Burr distribution")
    start <- startarg_fellerpareto_family(1/x, "burr")
  }else if(distr == "pareto4")
  {
    muhat <- murel(x)
    alphastar <- 2
    q1 <- quantile(x, probs=1/4, na.rm=TRUE)
    q3 <- quantile(x, probs=3/4, na.rm=TRUE)
    gammahat <- as.numeric(gammarelP4(alphastar, q1, q3))
    thetahat <- as.numeric(thetarelP4(alphastar, gammahat, q1, q3))
    alphahat <- alpharelP4(x, muhat, thetahat, gammahat)
    
    start <- list(min=muhat, shape1=alphahat, shape2=gammahat, scale=thetahat)
  }else if(distr == "pareto3")
  {
    muhat <- murel(x)
    alphastar <- 1 #true value in that case
    q1 <- quantile(x, probs=1/4, na.rm=TRUE)
    q3 <- quantile(x, probs=3/4, na.rm=TRUE)
    gammahat <- as.numeric(gammarelP4(alphastar, q1, q3))
    thetahat <- as.numeric(thetarelP4(alphastar, gammahat, q1, q3))
    
    start <- list(min=muhat, shape=gammahat, scale=thetahat)
    
  }else if(distr == "pareto2")
  {
    muhat <- murel(x)
    gammastar <- 1 #true value
    alphastar <- 2
    q1 <- quantile(x, probs=1/4, na.rm=TRUE)
    q3 <- quantile(x, probs=3/4, na.rm=TRUE)
    thetahat <- as.numeric(thetarelP4(alphastar, gammastar, q1, q3))
    alphahat <- alpharelP4(x, muhat, thetahat, gammastar)
    
    start <- list(min=muhat, shape=alphahat, scale=thetahat)
  }else if (distr == "pareto1")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a Pareto distribution")
    muhat <- murel(x)
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
    thetahat <- as.numeric(thetarelP4(alphastar, gammastar, q1, q3))
    alphahat <- alpharelP4(x, mustar, thetahat, gammastar)
    
    start <- list(shape=alphahat, scale=thetahat)
  }else if (distr == "invpareto")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a inverse Pareto  distribution")
    start <- startarg_fellerpareto_family(1/x, "pareto")
  }else if (distr == "llogis")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-logistic  distribution")
    q25 <- as.numeric(quantile(x, 0.25, na.rm=TRUE))
    q75 <- as.numeric(quantile(x, 0.75, na.rm=TRUE))
    shape <- 2*log(2)/(log(q75)-log(q25))
    scale <- exp(log(q75)+log(q25))/2
    start <- list(shape=shape, scale=scale)
  }else if (distr == "paralogis")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a paralogistic  distribution")
    alphahat <- alpharelPL(x)
    thetahat <- thetarelPL(x, alphahat)
    start <- list(shape = alphahat, scale = thetahat)
  }else if (distr == "invparalogis")
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a inverse paralogistic  distribution")
    start <- startarg_fellerpareto_family(1/x, "paralogis")
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
    if(v > 0)
    {
      alphahat <- m^2/v
      thetahat <- v/m
    }else #exponential case
    {
      alphahat <- 1 
      thetahat <- m
    } 
    start <- list(shape=alphahat, rate=1/thetahat)
  }else if (distr == "weibull") 
  {
    if (any(x < 0)) 
      stop("values must be positive to fit a Weibull distribution")
    #q25 <- as.numeric(quantile(x, 0.25, na.rm = TRUE))
    #q75 <- as.numeric(quantile(x, 0.75, na.rm = TRUE))
    #if(q25 < q75) #check to avoid division by zero
    #{
    #  q50 <- median(x, na.rm = TRUE)
    #  tauhat <- (log(-log(1-1/4)) - log(-log(1-3/4))) / (log(q25) - log(q75))
    #  thetahat <- q50/(-log(1-1/2))^(tauhat)
    #}else
    #{
      m <- mean(log(x), na.rm = TRUE)
      v <- var(log(x), na.rm = TRUE)
      tauhat <- 1.2/sqrt(v)
      thetahat <- exp(m + 0.572/tauhat)
    #}
    
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