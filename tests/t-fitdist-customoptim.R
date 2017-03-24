myoptimize <- function(fn,par,ui,ci,...){
  res <- constrOptim(f=fn,theta=par,method="Nelder-Mead",ui=ui,ci=ci, ...)
  standardres <- c(res,convergence=0,value=res$objective,par=res$minimum,hessian=NA)
  return(standardres)
}


#one parameter example
x <- rexp(1000)
library(fitdistrplus)
fitdist(x, "exp")

#binding example 
fitdist(x, "exp", custom.optim=myoptimize, ui=1, ci=2, start=list(rate=10))
fitdist(x, "exp", lower= 2, optim.method="L-BFGS-B")

#two parameter example
x <- rbeta(100, pi, 1/pi)

fitdist(x, "beta")

#binding example 
fitdist(x, "beta", custom.optim=myoptimize, ui=rbind(1,1), ci=c(1/2,1/2), start=list(shape1=5, shape2=5))
fitdist(x, "beta", lower= c(1/2,1/2), optim.method="L-BFGS-B")



#true example
library(GeneralizedHyperbolic)
args(dnig)


x <- rnig(1000, 3, 1/2, 1/2, 1/4)
ui<-rbind(c(0,1,0,0),c(0,0,1,0),c(0,0,1,-1),c(0,0,1,1))
ci<-c(0,0,0,0)
fitdist(x, "nig", custom.optim=myoptimize, ui=ui, ci=ci, start=list(mu = 0, delta = 1, alpha = 1, beta = 0))
?dnig
