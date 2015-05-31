library(fitdistrplus)
# (1) Fit of a lognormal distribution to bacterial contamination data
#
data(smokedfish)
fitsf  <-  fitdistcens(smokedfish,"lnorm")
summary(fitsf)
# default plot using the Turnbull algorithm
plot(fitsf)
# plot using the Turnbull algorithm with confidence intervals for the 
# empirical distribution
plot(fitsf, Turnbull.confint = TRUE)
# basic plot using intervals and points (see ?plotdiscens for details)
plot(fitsf, Turnbull = FALSE)
# plot of the same fit using the Turnbull algorithm in logscale
cdfcompcens(fitsf,main="bacterial contamination fits",
            xlab="bacterial concentration (CFU/g)",ylab="F",
            addlegend = FALSE,lines01 = TRUE, xlogscale = TRUE, xlim = c(1e-2,1e2))
# zoom on large values of F
cdfcompcens(fitsf,main="bacterial contamination fits",
            xlab="bacterial concentration (CFU/g)",ylab="F",
            addlegend = FALSE,lines01 = TRUE, xlogscale = TRUE, 
            xlim = c(1e-2,1e2),ylim=c(0.4,1))

# (2) Fit of a normal distribution on acute toxicity values 
# of fluazinam (in decimal logarithm) for
# macroinvertebrates and zooplancton, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology
#

data(fluazinam)
log10EC50 <-log10(fluazinam)
fln <- fitdistcens(log10EC50,"norm")
fln
summary(fln)
plot(fln)

# (3) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view dedicated to 
# probability distributions
#

dgumbel  <-  function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel  <-  function(q,a,b) exp(-exp((a-q)/b))
qgumbel  <-  function(p,a,b) a-b*log(-log(p))
fg <- fitdistcens(log10EC50,"gumbel",start=list(a=1,b=1))
summary(fg)
plot(fg)

# (4) comparison of fits of various distributions
# 

fll <- fitdistcens(log10EC50,"logis")
summary(fll)

cdfcompcens(list(fln,fll,fg),legendtext=c("normal","logistic","gumbel"),
            xlab = "log10(EC50)")

# (5) how to change the optimisation method?
#

fitdistcens(log10EC50,"logis",optim.method="Nelder-Mead")
fitdistcens(log10EC50,"logis",optim.method="BFGS") 
fitdistcens(log10EC50,"logis",optim.method="SANN") 

# (6) custom optimisation function - example with the genetic algorithm
#
  #wrap genoud function rgenoud package
  mygenoud  <-  function(fn, par, ...) 
  {
    require(rgenoud)
    res  <-  genoud(fn, starting.values=par, ...)        
    standardres  <-  c(res, convergence=0)
    
    return(standardres)
  }
  
  # call fitdistcens with a 'custom' optimization function
  fit.with.genoud <- fitdistcens(log10EC50, "logis", custom.optim=mygenoud, nvars=2,    
                                 Domains=cbind(c(0,0), c(5, 5)), boundary.enforcement=1, 
                                 print.level=1, hessian=TRUE)
  
  summary(fit.with.genoud)


# (7) estimation of the mean of a normal distribution 
# by maximum likelihood with the standard deviation fixed at 1 using the argument fix.arg
#
flnb <- fitdistcens(log10EC50, "norm", start = list(mean = 1),fix.arg = list(sd = 1))

# (8) Fit of a lognormal distribution on acute toxicity values of fluazinam for
# macroinvertebrates and zooplancton, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5 percent quantile value of 
# the fitted distribution (which is called the 5 percent hazardous concentration, HC5,
# in ecotoxicology) and estimation of other quantiles.

data(fluazinam)
log10EC50 <-log10(fluazinam)
fln <- fitdistcens(log10EC50,"norm")

quantile(fln, probs = 0.05)
quantile(fln, probs = c(0.05, 0.1, 0.2))


# (9) check keepdata
#
 x <- rexp(1e3, 5)
# x <- data.frame(left=x, right=x+rexp(x, 1/2))
x <- data.frame(left=x, right=x)
f1 <- fitdistcens(x, "exp", keepdata=FALSE)
f2 <- fitdistcens(x, "exp", keepdata=TRUE)
 
plot(f1)
plot(f2)  
