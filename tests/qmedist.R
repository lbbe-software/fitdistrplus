library(fitdistrplus)

# (1) basic fit of a normal distribution 
#

set.seed(1234)
x1 <- rnorm(n=100)
qmedist(x1, "norm", probs=c(1/3, 2/3))


# (2) defining your own distribution functions, here for the Gumbel 
# distribution for other distributions, see the CRAN task view dedicated 
# to probability distributions

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
qgumbel <- function(p, a, b) a - b*log(-log(p))
qmedist(x1, "gumbel", probs=c(1/3, 2/3), start=list(a=10,b=5))

# (3) fit a discrete distribution (Poisson)
#

set.seed(1234)
x2 <- rpois(n=30,lambda = 2)
qmedist(x2, "pois", probs=1/2)

# (4) fit a finite-support distribution (beta)
#

set.seed(1234)
x3 <- rbeta(n=100,shape1=5, shape2=10)
qmedist(x3, "beta", probs=c(1/3, 2/3))


# (5) fit frequency distributions on USArrests dataset.
#

x4 <- USArrests$Assault
qmedist(x4, "pois", probs=1/2)
qmedist(x4, "nbinom", probs=c(1/3, 2/3))

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
x <- c(rnorm(1000, 5),  rnorm(1000, 10))
#QME
fit2 <- qmedist(x, "norm2", probs=c(1/6, 1/4, 1/3, 1/2, 2/3), 
	start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2), 
	lower=c(0, 0, 0, 0, 0), upper=c(1/2, Inf, Inf, Inf, Inf))



# (7) test error messages
#

dnorm3 <- qnorm3 <- function(x, a)
  "NA"
x <- rexp(10)

#should get a one-line error 
res <- qmedist(x, "norm3", start=list(a=1), probs=1/2)
#as in 
attr(try(log("a"), silent=TRUE), "condition")


# (8) weighted QME
#
n <- 1e6
n <- 1e2
x <- rpois(n, 10)
xtab <- table(x)
xval <- sort(unique(x))
f1 <- qmedist(x, "pois", start=list(lambda=mean(x)), lower=0, upper=100, probs=1/2) #, control=list(trace=1, REPORT=1)
f2 <- qmedist(xval, "pois", weights=xtab, start=list(lambda=mean(x)), probs=1/2)

f1$estimate
f2$estimate #should be identical

x <- rexp(n)
f3 <- qmedist(x, "exp", probs=1/2)
f4 <- qmedist(x, "exp", weights=c(rep(1, n/2), sqrt(1:(n/2))), probs=1/2)
f3$estimate
f4$estimate

f3$loglik
f4$loglik

median(x)
median(tail(x, 50))

# (9) test the component optim.message
x <- rnorm(1000)
#change parameter to obtain unsuccessful convergence
qmedist(x, "norm", probs=1:2/3, control=list(maxit=2), start=list(mean=1e5, sd=1), optim.method="L-BFGS-B", lower=0)

