library(fitdistrplus)




#(1) beta distribution
#

n <- 100
set.seed(12345)
x <- rbeta(n, 3, 3/4)

psi <- function(x) digamma(x)
grbetalnl <- function(x, a, b)
  c(log(x)-psi(a)+psi(a+b), log(1-x)-psi(b)+psi(a+b))
#grbetalnl(x, 3, 4)
grlnL <- function(par, obs, ...) 
  -rowSums(sapply(obs, function(x) grbetalnl(x, a=par[1], b=par[2])))
#rowSums(sapply(x, function(x) grbetalnl(x, 3, 4)))
#grlnL(c(3, 4), x)
#grlnL(c(3, 3/4), x)


constrOptim2 <- function(par, fn, gr=NULL, ui, ci, ...)
  constrOptim(theta=unlist(par), f=fn, grad=gr, ui=ui, ci=ci, ...)



#control parameters
ctr <- list(trace=3, REPORT=1, maxit=1000)
ctr <- list(trace=0, REPORT=1, maxit=1000)

bfgs_gr$time <- system.time(bfgs_gr <- mledist(x, dist="beta", optim.method="BFGS", gr=grlnL, control=ctr))[3]
bfgs <- mledist(x, dist="beta", optim.method="BFGS", control=ctr)

lbfgs_gr <- mledist(x, dist="beta", optim.method="L-BFGS-B", gr=grlnL, control=ctr, lower=c(0,0))
lbfgs <- mledist(x, dist="beta", optim.method="L-BFGS-B", control=ctr, lower=c(0,0))

cg_gr <- mledist(x, dist="beta", optim.method="CG", gr=grlnL, control=ctr)
cg <- mledist(x, dist="beta", optim.method="CG", control=ctr)

nm_gr <- mledist(x, dist="beta", optim.method="Nelder", gr=grlnL, control=ctr)
nm <- mledist(x, dist="beta", optim.method="Nelder", control=ctr)

constr_nm_gr <- mledist(x, dist="beta", custom.optim=constrOptim2,
      ui = diag(2), ci = c(0, 0), optim.method="Nelder", gr=grlnL, control=ctr)
constr_nm <- mledist(x, dist="beta", custom.optim=constrOptim2,
      ui = diag(2), ci = c(0, 0), optim.method="Nelder", control=ctr)

constr_bfgs_gr <- mledist(x, dist="beta", custom.optim=constrOptim2,
        ui = diag(2), ci = c(0, 0), optim.method="BFGS", gr=grlnL, control=ctr)
constr_bfgs <- mledist(x, dist="beta", custom.optim=constrOptim2,
        ui = diag(2), ci = c(0, 0), optim.method="BFGS", control=ctr)

constr_cg_gr <- mledist(x, dist="beta", custom.optim=constrOptim2,
                        ui = diag(2), ci = c(0, 0), optim.method="CG", gr=grlnL, control=ctr)
constr_cg <- mledist(x, dist="beta", custom.optim=constrOptim2,
                     ui = diag(2), ci = c(0, 0), optim.method="CG", control=ctr)

lnL <- function(par, fix.arg, obs, ddistnam) 
{
  fitdistrplus:::loglikelihood(par, fix.arg, obs, ddistnam) 
}
  

constrOptim2(c(shape1=1, shape2=1), lnL, obs=x, fix.arg=NULL, ddistnam="dbeta",
             ui = diag(2), ci = c(0, 0))


#no log
dbeta3 <- function(x, shape1, shape2)
  dbeta(x, shape1, shape2)

#Ripley trick : param transform
dbeta2 <- function(x, shape1, shape2, log)
  dbeta(x, exp(shape1), exp(shape2), log=log)
pbeta2 <- function(q, shape1, shape2, log.p)
  pbeta(q, exp(shape1), exp(shape2), log.p=log.p)

bfgs2 <- mledist(x, dist="beta2", optim.method="BFGS", control=ctr, 
                start=list(shape1=0, shape2=0))
bfgs3 <- mledist(x, dist="beta3", optim.method="BFGS", control=ctr, 
                 start=list(shape1=1, shape2=1))




getval <- function(x)
  c(x$estimate, loglik=x$loglik, x$counts)
getval2 <- function(x)
  c(exp(x$estimate), loglik=x$loglik, x$counts)

cbind(trueval=c(3, 3/4, lnL(c(3, 3/4), NULL, x, "dbeta"), NA, NA),
      NM=getval(nm), NMgrad=getval(nm_gr), 
      constr_NM=getval(constr_nm), constr_NMgrad=getval(constr_nm_gr),
      CG=getval(cg), CGgrad=getval(cg_gr), 
      constr_CG=getval(constr_cg), constr_CGgrad=getval(constr_cg_gr),
      BFGS=getval(bfgs), BFGSgrad=getval(bfgs_gr),
      constr_BFGS=getval(constr_bfgs), constr_BFGSgrad=getval(constr_bfgs_gr),
      BFGS_exp=getval2(bfgs2), BFGS_nolog=getval(bfgs3))



llsurface(min.arg=c(0.1, 0.1), max.arg=c(7, 3), plot.arg=c("shape1", "shape2"),
          lseq=50, data=x, distr="beta")
points(bfgs$estimate[1], bfgs$estimate[2], pch="+", col="red")
points(3, 3/4, pch="x", col="green")


