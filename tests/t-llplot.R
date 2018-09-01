require(fitdistrplus)

# (1) tests with the Burr distribution (three parameters)
#
if(any(installed.packages()[, "Package"] == "actuar"))
{
  require(actuar)
  
  data(endosulfan)
  ATV <-endosulfan$ATV
  library("actuar")
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
n <- 1e6
n <- 1e2
x <- rpois(n, 10)
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
if(FALSE)
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