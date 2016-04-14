require(fitdistrplus)

x <- rgamma(1e3, shape=3/2, rate= 1/2)

initval <- unlist(fitdistrplus:::start.arg.default(x, "gamma"))

fitdistrplus:::fitbench(x, "gamma", "mle")

fitdistrplus:::fitbench(x, "gamma", "mle", lower=0)

#compute gradient of the log-likelihood
psi <- function(x) digamma(x)
#individual contribution
grgamlnl <- function(x, shape, rate) 
  c(log(x)-log(rate)-psi(shape), x/rate^2-shape/rate)
#total grad loglik
grgam <- function(par, obs, ...) 
{
  n <- length(obs)
  res <- grgamlnl(obs, shape=par[1], rate=par[2])
  c(-sum(res[1:n]), -sum(res[1:n+n]))
}  
  

# grgam(1:2, x)
# grgamlnl(x[1], shape=1, rate=2)
# grgamlnl(x[1:2], shape=1, rate=2)

fitdistrplus:::fitbench(x, "gamma", "mle", grad=grgam)

fitdistrplus:::fitbench(x, "gamma", "mle", grad=grgam, lower=0)

#mledist(x, "gamma", grad=grgam, lower=0, optim.method = "CG")
#mledist(x, "gamma", grad=grgam, lower=0, optim.method = "BFGS")
