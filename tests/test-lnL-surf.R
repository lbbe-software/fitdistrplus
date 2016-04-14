library(fitdistrplus)



#(1) beta distribution
#

n <- 100
x <- rbeta(n, 3, 3/4)

llsurface(data = x, distr = "beta", plot.arg=c("shape1", "shape2"),
          min.arg=c(0.1, 0.1), max.arg=c(7, 3))
llsurface(data = x, distr = "beta", plot.arg=c("shape1", "shape2"),
          min.arg=c(0.1, 0.1), max.arg=c(7, 3), col = FALSE )
points(3, 3/4, pch="+", col="red")


llcurve(data = x, distr = "beta", plot.arg = "shape1", min.arg = 0.1, max.arg = 7, 
        fix.arg = list(shape2 = 3/4), lseq=100, col = "blue")
llcurve(data = x, distr = "beta", plot.arg = "shape2", min.arg = 0.1, max.arg = 7, 
        fix.arg = list(shape1 = 3), lseq=100, col = "red")




#test

psi <- function(x) digamma(x)

grbetalnl <- function(x, a, b)
  c(log(x)-psi(a)+psi(a+b), log(1-x)-psi(b)+psi(a+b))

grbetalnl(x, 3, 4)

grlnL <- function(par, obs, ...) 
  -rowSums(sapply(obs, function(x) grbetalnl(x, a=par[1], b=par[2])))

rowSums(sapply(x, function(x) grbetalnl(x, 3, 4)))
grlnL(c(3, 4), x)
grlnL(c(3, 3/4), x)


ctr <- list(trace=0, REPORT=1, maxit=1000)

bfgs_gr <- mledist(x, dist="beta", optim.method="BFGS", gr=grlnL, control=ctr)
bfgs <- mledist(x, dist="beta", optim.method="BFGS", control=ctr)

cg_gr <- mledist(x, dist="beta", optim.method="CG", gr=grlnL, control=ctr)
cg <- mledist(x, dist="beta", optim.method="CG", control=ctr)

nm_gr <- mledist(x, dist="beta", optim.method="Nelder", gr=grlnL, control=ctr)
nm <- mledist(x, dist="beta", optim.method="Nelder", control=ctr)

getval <- function(x)
  c(x$estimate, loglik=x$loglik, x$counts)

cbind(NM=getval(nm), NMgrad=getval(nm_gr), CG=getval(cg), 
      CGgrad=getval(cg_gr), BFGS=getval(bfgs), BFGSgrad=getval(bfgs_gr))



llsurface(data = x, distr = "beta", plot.arg = c("shape1", "shape2"),
          min.arg = c(0.1, 0.1), max.arg = c(7, 3), col.pal = heat.colors(50))
points(bfgs$estimate[1], bfgs$estimate[2], pch="+", col="red")
points(3, 3/4, pch="x", col="green")


