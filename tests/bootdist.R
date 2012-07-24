library(fitdistrplus)

# (1) basic fit of a normal distribution with maximum likelihood estimation
# followed by parametric bootstrap
#
x1 <- c(6.4, 13.3, 4.1, 1.3, 14.1, 10.6, 9.9, 9.6, 15.3, 22.1, 13.4, 
13.2, 8.4, 6.3, 8.9, 5.2, 10.9, 14.4)
f1 <- fitdist(x1, "norm", method="mle")
b1 <- bootdist(f1)
print(b1)
plot(b1)
summary(b1)

# (2) non parametric bootstrap
#
b1np <- bootdist(f1, bootmethod="nonparam", niter=101)
summary(b1np)

# (3) fit of a gamma distribution followed by parametric bootstrap
#
f1b <- fitdist(x1, "gamma", method="mle")
b1b <- bootdist(f1b, niter=101)
summary(b1b)

# (4) fit of a gamma distribution with control of the optimization
# method,  followed by parametric bootstrap
#
f1c <- fitdist(x1, "gamma", optim.method="L-BFGS-B", lower=c(0, 0))
b1c <- bootdist(f1c, niter=101)
summary(b1c)

# (5) estimation of the standard deviation of a normal distribution 
# by maximum likelihood with the mean fixed at 10 using the argument fix.arg
# followed by parametric bootstrap
#
f1d  <- fitdist(x1, "norm", start=list(sd=5), fix.arg=list(mean=10))
b1d <- bootdist(f1d, niter=101)
summary(b1d)
plot(b1d)

# (6) fit of a discrete distribution by matching moment estimation 
# (using a closed formula) followed by parametric bootstrap
#
x2 <- c(rep(4, 1), rep(2, 3), rep(1, 7), rep(0, 12))
f2 <- fitdist(x2, "pois", method="mme")
b2 <- bootdist(f2, niter=101)
plot(b2, pch=16)
summary(b2)

# (7) fit of a Weibull distribution to serving size data by maximum likelihood 
# estimation or by quantile matching estimation (in this example matching 
# first and third quartiles) followed by parametric bootstrap
#
data(groundbeef)
serving <- groundbeef$serving

fWmle <- fitdist(serving, "weibull")
bWmle <- bootdist(fWmle, niter=101)
summary(bWmle)

fWqme <- fitdist(serving, "weibull", method="qme", probs=c(0.25, 0.75))
bWqme <- bootdist(fWqme, niter=101)
summary(bWqme)

# (8) Fit of a Pareto distribution by numerical moment matching estimation
# followed by parametric bootstrap
#
if(any(installed.packages()[, "Package"] == "actuar"))
{
    require(actuar)
#simulate a sample
    x4 <- rpareto(1000,  6,  2)
    memp <- function(x,  order)
	ifelse(order == 1,  mean(x),  sum(x^order)/length(x))
	
    f4 <- fitdist(x4,  "pareto",  "mme",  order=1:2,  
				  start=c(shape=10,  scale=10),  
				  lower=1,  memp="memp",  upper=50)
	
    b4 <- bootdist(f4,  niter=101)
    summary(b4)
	
    b4npar <- bootdist(f4,  niter=101,  bootmethod="nonparam")
    summary(b4npar)
}	


# (9) Fit of a uniform distribution using Cramer-von Mises 
# followed by parametric boostrap
# 

u <- runif(50, min=5, max=10)

fu <- fitdist(u, "unif", method="mge", gof="CvM")
bu <- bootdist(fu,  bootmethod="param", niter=101)
summary(bu)
plot(bu)

