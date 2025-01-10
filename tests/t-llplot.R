require("fitdistrplus")
visualize <- FALSE # TRUE for manual tests with visualization of results
nsample <- 10000
nsample <- 10

# (1) tests with the Burr distribution (three parameters)
#
if(any(installed.packages()[, "Package"] == "actuar"))
{
  require(actuar)
  
  data(endosulfan)
  ATV <-endosulfan$ATV
  require("actuar")
  fBurr <- fitdist(ATV, "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1))
  llplot(fBurr)
  
  fBurr2 <- fitdist(ATV, "burr", start = list(shape1 = 0.3, shape2 = 1),
                    fix.arg = list(rate = 1.5))
  llplot(fBurr2)
  
  fBurr3 <- fitdist(ATV, "burr", start = list(shape1 = 0.3, rate = 1),
                    fix.arg = list(shape2 = 1.5))
  llplot(fBurr3)
}


# (2) An example on discrete data with or without weights
#
set.seed(1234)
x <- rpois(nsample, 10)
xtab <- table(x)
xval <- sort(unique(x))
f1 <- fitdist(x, "pois")
f2 <- fitdist(xval, "pois", weights = xtab)

f1$estimate
f2$estimate # should give the same
llplot(f1, fit.show = TRUE)
llplot(f2, fit.show = TRUE) # should give the same
llplot(f1, loglik = FALSE, fit.show = TRUE)
llplot(f2, loglik = FALSE,fit.show = TRUE) # should give the same

# (3) An example on censored data with or without weights
#
if(visualize)
{
  data(salinity)
  salinity.unique <- unique(salinity)
  string.unique <- paste(salinity.unique$left, salinity.unique$right)
  string.salinity <- paste(salinity$left, salinity$right)
  nobs <- nrow(salinity.unique)
  salinity.weights <- numeric(nobs)
  for (i in 1:nobs)
  {
    salinity.weights[i] <- length(which(string.salinity == string.unique[i]))
  }
  cbind(salinity.unique, salinity.weights)
  
  (fa <- fitdistcens(salinity, "lnorm"))
  (fb <- fitdistcens(salinity.unique, "lnorm", weights = salinity.weights))
  llplot(fa, fit.show = TRUE)
  llplot(fb, fit.show = TRUE) # should give the same
  llplot(fa, fit.show = TRUE, loglik = FALSE)
  llplot(fb, fit.show = TRUE, loglik = FALSE) # should give the same
  
}

# (4) An example with NaN stderror
# 

if(visualize)
{
  claims <- read.csv("~/Documents/recherche-enseignement/code/R/riskassessment/bug/20250107/Claims.csv")
  x <- claims$UltimateCost/1000
  
  fit_B_mle <- fitdist(x, "burr", method="mle", lower=0)
  llplot(fit_B_mle, expansion = 10, fit.show = TRUE)
  
  fit_IB_mle <- fitdist(x, "invburr", method="mle", lower=0) # converges to a wrong solution
  llplot(fit_IB_mle, expansion = 10, fit.show = TRUE)
  
  fit_IB_mle$estimate[1] <- NaN
  llplot(fit_IB_mle, expansion = 10, fit.show = TRUE)
}