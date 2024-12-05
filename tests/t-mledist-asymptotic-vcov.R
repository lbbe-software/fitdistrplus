require("fitdistrplus")
nsample <- 1e3
nsample <- 10

#### (1) Gamma example ####

truetheta <- c("alpha"=0.19, "beta"=5.18)
x <- rgamma(nsample, truetheta["alpha"], truetheta["beta"])
f1 <- mledist(x, "gamma", calcvcov = TRUE)
f1$vcov
f1$estimate
infoFisher <- function(alpha, beta)
{
  cbind(c(trigamma(alpha), -1/beta),
        c(-1/beta, alpha/beta))
}
solve(infoFisher(0.19, 5.18))/nsample


if(FALSE)
{
  #check with MASS::fitdistr()
  
  mledist(rgamma(1e2, truetheta["alpha"], truetheta["beta"]), "gamma", calcvcov = TRUE)$vcov
  mledist(rgamma(1e3, truetheta["alpha"], truetheta["beta"]), "gamma", calcvcov = TRUE)$vcov
  mledist(rgamma(1e4, truetheta["alpha"], truetheta["beta"]), "gamma", calcvcov = TRUE)$vcov
  
  
  MASS::fitdistr(rgamma(1e2, truetheta["alpha"], truetheta["beta"]), "gamma")$vcov
  MASS::fitdistr(rgamma(1e3, truetheta["alpha"], truetheta["beta"]), "gamma")$vcov
  MASS::fitdistr(rgamma(1e4, truetheta["alpha"], truetheta["beta"]), "gamma")$vcov
  
  
}


# (2) fit a Pareto distribution
#

if(any(installed.packages()[, "Package"] == "actuar"))
{
  require("actuar")
  #simulate a sample
  x4 <- rpareto(nsample, 6, 2)
  
  #empirical raw moment
  memp <- function(x, order)
    mean(x^order)
  
  #fit
  res <- mledist(x4, "pareto", start=list(shape=10, scale=10), 
          lower=1, upper=Inf, calcvcov = TRUE)
}

# (3) truncated distribution
#

dtiexp <- function(x, rate, low, upp)
{
  PU <- pexp(upp, rate=rate, lower.tail = FALSE)
  PL <- pexp(low, rate=rate)
  dexp(x, rate) * (x >= low) * (x <= upp) + PL * (x == low) + PU * (x == upp)
}
ptiexp <- function(q, rate, low, upp)
  pexp(q, rate) * (q >= low) * (q <= upp) + 1 * (q > upp)
n <- 100; x <- pmax(pmin(rexp(n), 3), .5)

mledist(x, "tiexp", start=list(rate=3, low=0, upp=20), calcvcov = TRUE)
