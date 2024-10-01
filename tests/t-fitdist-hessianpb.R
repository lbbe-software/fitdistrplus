require("fitdistrplus")



#an example with lower bounds for gamma

data(groundbeef)
serving <- groundbeef$serving
#(f1 <-  fitdist(serving, "gamma", optim.method="Nelder-Mead"))
(f1 <-  fitdist(serving, "gamma", optim.method="Nelder-Mead", lower=0))

f1 <- mledist(serving, "gamma", optim.method="Nelder-Mead", lower=0)

f1 <- fitdist(serving, "gamma", start=as.list(f1$estimate))

#manual optimization 
lnl <- function(theta, obs, dfun)
  sum(do.call(dfun, c(list("x"=obs), theta, "log"=TRUE)))

H <- optimHess(coef(f1), lnl, obs=serving, dfun=dgamma)

-solve(H)
vcov(fitdist(serving, "gamma", start=as.list(f1$estimate)))



#an example with lower bounds for ztbinom
if(FALSE)
{
    require("actuar")
  n <- 1e3
  x <- rztbinom(n, 30, 1/2)
  lnl(c(35, 1/2), obs=x, dfun=dztbinom)
  
  optim(c(35, 1/2), lnl, obs=x, dfun=dztbinom, control=list(fnscale=-1))
  optimHess(c(35, ), lnl, obs=x, dfun=dztbinom)
  
  mledist(x, "ztbinom", fix.arg = list(size=30))
  fitdist(x, "ztbinom", fix.arg = list(size=30))
  
  mledist(x, "ztbinom", lower=0)
  fitdist(x, "ztbinom", lower=0)
  
  dztbinom(5, 30.5, 1/2)
  
  llsurface(data = x, distr="ztbinom", plot.arg=c("size", "prob"), 
            min.arg=c(max(x)+1, 0.4), max.arg=c(31, 0.6), back.col = FALSE,
            main="log-likelihood for ZT binomial distribution")
  
}