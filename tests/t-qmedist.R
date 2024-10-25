require("fitdistrplus")
nsample <- 10

# (1) basic fit of a normal distribution 
#

set.seed(1234)
x1 <- rnorm(n=nsample)
qmedist(x1, "norm", probs=c(1/3, 2/3))

#fitted coef is -0.1486435  0.9892013, fitted loglik is -13.92195

# (2) defining your own distribution functions, here for the Gumbel 
# distribution for other distributions, see the CRAN task view dedicated 
# to probability distributions

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
qgumbel <- function(p, a, b) a - b*log(-log(p))
qmedist(x1, "gumbel", probs=c(1/3, 2/3), start=list(a=10,b=5))

#fitted coef is -0.4935571  0.8557281, fitted loglik is -16.79828

# (3) fit a discrete distribution (Poisson)
#

set.seed(1234)
x2 <- rpois(n=nsample,lambda = 2)
qmedist(x2, "pois", probs=1/2)

#fitted coef is 1.7, fitted loglik is -15.31626

# (4) fit a finite-support distribution (beta)
#

set.seed(1234)
x3 <- rbeta(n=nsample, shape1=5, shape2=10)
qmedist(x3, "beta", probs=c(1/3, 2/3))

#fitted coef is 4.010999 8.397853, fitted loglik is 6.642124

# (5) fit frequency distributions on USArrests dataset.
#

x4 <- USArrests$Assault
qmedist(x4, "pois", probs=1/2)
qmedist(x4, "nbinom", probs=c(1/3, 2/3))

#fitted coef is 2.518966 182.313344, fitted loglik is 292.5969

# (6) normal mixture
#

#mixture of two normal distributions
#density
dnorm2 <- function(x, poid, m1, s1, m2, s2)
	poid*dnorm(x, m1, s1) + (1-poid)*dnorm(x, m2, s2)
#numerically approximate quantile function
qnorm2 <- function(p, poid, m1, s1, m2, s2)
{
	L2 <- function(x, prob)
		(prob - pnorm2(x, poid, m1, s1, m2, s2))^2	
	sapply(p, function(pr) optimize(L2, c(-20, 30), prob=pr)$minimum)
}	
#distribution function		
pnorm2 <- function(q, poid, m1, s1, m2, s2)
	poid*pnorm(q, m1, s1) + (1-poid)*pnorm(q, m2, s2)		


#basic normal distribution
x <- c(rnorm(nsample, 5),  rnorm(nsample, 10))
#QME
fit2 <- qmedist(x, "norm2", probs=c(1/6, 1/4, 1/3, 1/2, 2/3), 
	start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2), 
	lower=c(0, 0, 0, 0, 0), upper=c(1/2, Inf, Inf, Inf, Inf))

fit2
#fitted coef is 0.3433528 4.2872449 0.3891135 9.2598612 1.7948554, fitted loglik is -38.14106

# (7) test error messages
#

dnorm3 <- qnorm3 <- function(x, a)
  "NA"
x <- rexp(nsample)

#should get a one-line error 
res <- qmedist(x, "norm3", start=list(a=1), probs=1/2)
#as in 
attr(try(log("a"), silent=TRUE), "condition")





try(qmedist(c(x1, "a"), "gamma", probs=1/2))
try(qmedist(c(x1, NA), "gamma", probs=1/2))
try(qmedist(c(x1, Inf), "gamma", probs=1/2))
try(qmedist(c(x1, -Inf), "gamma", probs=1/2))
try(qmedist(c(x1, NaN), "gamma", probs=1/2))


# (8) weighted QME
#
n <- 1e6
x <- rpois(nsample, 10)
xtab <- table(x)
xval <- sort(unique(x))
f1 <- qmedist(x, "pois", start=list(lambda=mean(x)), lower=0, upper=100, probs=1/2) #, control=list(trace=1, REPORT=1)
f2 <- qmedist(xval, "pois", weights=xtab, start=list(lambda=mean(x)), probs=1/2)

f1$estimate
f2$estimate #should be identical

x <- rexp(nsample)
f3 <- qmedist(x, "exp", probs=1/2)
w4 <- c(rep(1, nsample/2), round(sqrt(1:(nsample/2))))
f4 <- qmedist(x, "exp", weights=w4, probs=1/2)
c(f3$estimate, f4$estimate)

c(f3$loglik, f4$loglik)

#fitted coef is 0.4816191, fitted loglik is -16.95355


#try non integer weights
try(qmedist(x, "exp", weights=c(rep(1, n/2), sqrt(1:(n/2))), probs=1/2))

# (9) test the component optim.message
x <- rnorm(nsample)
#change parameter to obtain unsuccessful convergence
qmedist(x, "norm", probs=1:2/3, control=list(maxit=2), start=list(mean=1e5, sd=1), optim.method="L-BFGS-B", lower=0)

# (10) test bounds
x <- rnorm(nsample)
qmedist(x, "norm", probs=1:2/3, optim.method="L-BFGS-B", lower=c(-Inf, 0)) #via optim
qmedist(x, "norm", probs=1:2/3, optim.method="Nelder", lower=c(-Inf, 0)) #via constrOptim


# (11) large sample size issue

if(FALSE)
{
  set.seed(123)
  obs <- rlnorm(1e7, 3, 2)
  for(i in 2:7)
    cat(i, try(qmedist(obs[1:10^i], "lnorm", probs=1:2/3, control=list(trace=0, REPORT=1))$estimate, silent=TRUE), "\n")
  
  # 2 3.109257 1.767013 
  # 3 3.022615 1.857885 
  # 4 2.999376 1.978701 
  # 5 3.003344 2.00941 
  # 6 2.99881 2.001733 
  # 7 2.999859 2.000227 
}

