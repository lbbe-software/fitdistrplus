require("fitdistrplus")

nsample <- 500
nsample <- 10
visualize <- FALSE # TRUE for manual tests with visualization of results
set.seed(1234)


# (1) check input data
#
data(smokedfish)
fitsf  <-  fitdistcens(smokedfish, "lnorm")
summary(fitsf)

#fitted coef is -3.63    3.54, fitted loglik is -90.7

smokedfish2 <- tail(smokedfish, 20)

try(fitdistcens(rbind(smokedfish2, c(NA, NA)) , "lnorm"))

try(fitdistcens(rbind(smokedfish2, c(3, Inf)) , "lnorm"))

try(fitdistcens(rbind(smokedfish2, c(NaN, 3)) , "lnorm"))

try(fitdistcens(rbind(smokedfish2, c(5, 3)) , "lnorm"))

fitsf  <-  fitdistcens(smokedfish, "lnorm", calcvcov = FALSE)



# (6) custom optimisation function - example with the genetic algorithm
#
data(fluazinam)
log10EC50 <- log10(fluazinam)

summary(fitdistcens(log10EC50, "logis"))

#fitted coef is 2.152 0.691 , fitted loglik is -20.6

#wrap genoud function rgenoud package
mygenoud  <-  function(fn, par, ...) 
{
  require("rgenoud")
  res  <-  genoud(fn, starting.values=par, ...)        
  standardres  <-  c(res, convergence=0)
  
  return(standardres)
}

# call fitdistcens with a 'custom' optimization function
fit.with.genoud <- fitdistcens(log10EC50, "logis", custom.optim=mygenoud, nvars=2, 
                               start=list(location=1, scale=1),    
                               Domains=cbind(c(0,0), c(5, 5)), boundary.enforcement=1, 
                               print.level=1, hessian=TRUE)

summary(fit.with.genoud)


# (9) check keepdata
#
if (visualize) # LONG TO RUN ON CRAN AND NEEDS VISALIZATION OF RESULTS
{
  set.seed(1234)
  x <- rexp(1e3, 5)
  # x <- data.frame(left=x, right=x+rexp(x, 1/2))
  x <- data.frame(left=x, right=x)
  f1 <- fitdistcens(x, "exp", keepdata=FALSE)
  f2 <- fitdistcens(x, "exp", keepdata=TRUE)
  f1$censdata
  f2$censdata
  plot(f1)
  plot(f2)
  
  #fitted coef is 5, fitted loglik is 609
}


# (9) fixing parameters
#
x <- rexp(nsample, 5)
x <- data.frame(left=x, right=x+.1)

f1 <- fitdistcens(x, "gamma", fix.arg=list(shape=1.5))
f1
f1$fix.arg

#fitted coef is 1.5, fitted loglik is  -14.7

f1 <- fitdistcens(x, "gamma", fix.arg=function(x) list(shape=1.5))
f1
f1$fix.arg.fun


# (10) weights
#
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
(fb <- fitdistcens(salinity.unique, "lnorm", weights = salinity.weights)) # should give the same results


#fitted coef is 3.385   0.496, fitted loglik is  -139


# (11) check the warning messages when using weights in the fit followed by functions
# that do not yet take weights into account
# with an example to be used later to see if weights are well taken into account
#
x <- rexp(100, 5)
x <- sort(x)
x <- data.frame(left=x, right=x+.1)
(f <- fitdistcens(x, "gamma", weights=c(rep(10, 50), rep(1, 50))))
try(plot(f))
try(cdfcompcens(f))
(f2 <- fitdistcens(x, "weibull", weights=c(rep(10, 50), rep(1, 50))))
try(cdfcompcens(list(f, f2)))
try(bootdistcens(f))
