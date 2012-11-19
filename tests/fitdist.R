library(fitdistrplus)



# (1) basic fit of a gamma distribution by maximum likelihood estimation
#
data(groundbeef)
serving <- groundbeef$serving
fitg <- fitdist(serving, "gamma")
summary(fitg)
plot(fitg)
cdfcomp(fitg, addlegend=FALSE)


# (2) use the moment matching estimation (using a closed formula)
#

fitgmme <- fitdist(serving, "gamma", method="mme")
summary(fitgmme)
plot(fitgmme)
quantile(fitgmme)

# (3) fit and comparison of various fits
#
fitW <- fitdist(serving, "weibull")
fitg <- fitdist(serving, "gamma")
fitln <- fitdist(serving, "lnorm")
summary(fitW)
summary(fitg)
summary(fitln)
cdfcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))

# (4) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view 
# dedicated to probability distributions
#
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

fitgumbel <- fitdist(serving, "gumbel", start=list(a=10, b=10))
summary(fitgumbel)
plot(fitgumbel)

# (5) fit a discrete distribution (Poisson)
#

x2<-c(rep(4, 1), rep(2, 3), rep(1, 7), rep(0, 12))
f2<-fitdist(x2, "pois")
plot(f2)
summary(f2)
gofstat(f2)

# (6) how to change the optimisation method?
#

fitdist(serving, "gamma", optim.method="Nelder-Mead")
fitdist(serving, "gamma", optim.method="BFGS") 
fitdist(serving, "gamma", optim.method="SANN")

# (7) custom optimization function
#

#create the sample
mysample <- rexp(100, 5)
mystart <- list(rate=8)

res1 <- fitdist(mysample, dexp, start= mystart, optim.method="Nelder-Mead")

#show the result
summary(res1)

#the warning tell us to use optimise, because the Nelder-Mead is not adequate.

#to meet the standard 'fn' argument and specific name arguments, we wrap optimize, 
myoptimize <- function(fn, par, ...) 
{
    res <- optimize(f=fn, ..., maximum=FALSE)  
#assume the optimization function minimize
    
    standardres <- c(res, convergence=0, value=res$objective, 
					 par=res$minimum, hessian=NA)
    
    return(standardres)
}

#call fitdist with a 'custom' optimization function
res2 <- fitdist(mysample, dexp, start=mystart, custom.optim=myoptimize, 
interval=c(0, 100))

#show the result
summary(res2)


# (8) custom optimization function - another example with the genetic algorithm
#
if(any(installed.packages()[,"Package"] == "rgenoud"))
{

#set a sample
    x1 <- c(6.4, 13.3, 4.1, 1.3, 14.1, 10.6, 9.9, 9.6, 15.3, 22.1, 
			13.4, 13.2, 8.4, 6.3, 8.9, 5.2, 10.9, 14.4) 
    fit1 <- fitdist(x1, "gamma")
    summary(fit1)
	
#wrap genoud function rgenoud package
    mygenoud <- function(fn, par, ...) 
    {
        require(rgenoud)
        res <- genoud(fn, starting.values=par, ...)        
        standardres <- c(res, convergence=0)
		
        return(standardres)
    }
	
#call fitdist with a 'custom' optimization function
    fit2 <- fitdist(x1, "gamma", custom.optim=mygenoud, nvars=2,    
					Domains=cbind(c(0, 0), c(10, 10)), boundary.enforcement=1, 
					print.level=0, hessian=TRUE)
	
    summary(fit2)
}

# (9) estimation of the standard deviation of a normal distribution 
# by maximum likelihood with the mean fixed at 10 using the argument fix.arg
#
x1 <- c(6.4, 13.3, 4.1, 1.3, 14.1, 10.6, 9.9, 9.6, 15.3, 22.1, 
13.4, 13.2, 8.4, 6.3, 8.9, 5.2, 10.9, 14.4) 
fitdist(x1, "norm", start=list(sd=5), fix.arg=list(mean=10))

# (10) fit of a Weibull distribution to serving size data 
# by maximum likelihood estimation
# or by quantile matching estimation (in this example 
# matching first and third quartiles)
#
data(groundbeef)
serving <- groundbeef$serving

fWmle <- fitdist(serving, "weibull")
summary(fWmle)
plot(fWmle)
gofstat(fWmle)

fWqme <- fitdist(serving, "weibull", method="qme", probs=c(0.25, 0.75))
summary(fWqme)
plot(fWqme)
gofstat(fWqme)


# (11) Fit of a Pareto distribution by numerical moment matching estimation
#
if(any(installed.packages()[,"Package"] == "actuar"))
{

    require(actuar)
#simulate a sample
    x4 <- rpareto(1000, 6, 2)
	
#empirical raw moment
    memp <- function(x, order)
	ifelse(order == 1, mean(x), sum(x^order)/length(x))
	
	
#fit
    fP <- fitdist(x4, "pareto", method="mme", order=c(1, 2), memp="memp", 
				  start=c(shape=10, scale=10), lower=1, upper=Inf)
    summary(fP)
	
}

# (12) Fit of a Weibull distribution to serving size data by maximum 
# goodness-of-fit estimation using all the distances available
# 

data(groundbeef)
serving <- groundbeef$serving
(f1 <- fitdist(serving, "weibull", method="mge", gof="CvM"))
(f2 <- fitdist(serving, "weibull", method="mge", gof="KS"))
(f3 <- fitdist(serving, "weibull", method="mge", gof="AD"))
(f4 <- fitdist(serving, "weibull", method="mge", gof="ADR"))
(f5 <- fitdist(serving, "weibull", method="mge", gof="ADL"))
(f6 <- fitdist(serving, "weibull", method="mge", gof="AD2R"))
(f7 <- fitdist(serving, "weibull", method="mge", gof="AD2L"))
(f8 <- fitdist(serving, "weibull", method="mge", gof="AD2"))
cdfcomp(list(f1, f2, f3, f4, f5, f6, f7, f8))
cdfcomp(list(f1, f2, f3, f4, f5, f6, f7, f8), xlogscale=TRUE, xlim=c(8, 250), verticals=TRUE)

# (13) Fit of a uniform distribution using Cramer-von Mises or
# Kolmogorov-Smirnov distance
# 

u <- runif(50, min=5, max=10)

fuCvM <- fitdist(u, "unif", method="mge", gof="CvM")
summary(fuCvM)
plot(fuCvM)
gofstat(fuCvM)

fuKS <- fitdist(u, "unif", method="mge", gof="KS")
summary(fuKS)
plot(fuKS)
gofstat(fuKS)

# (14) scaling problem
#

x2 <- c(-0.00707717, -0.000947418, -0.00189753, 
-0.000474947, -0.00190205, -0.000476077, 0.00237812, 0.000949668, 
0.000474496, 0.00284226, -0.000473149, -0.000473373, 0, 0, 0.00283688, 
-0.0037843, -0.0047506, -0.00238379, -0.00286807, 0.000478583, 
0.000478354, -0.00143575, 0.00143575, 0.00238835, 0.0042847, 
0.00237248, -0.00142281, -0.00142484, 0, 0.00142484, 0.000948767, 
0.00378609, -0.000472478, 0.000472478, -0.0014181, 0, -0.000946522, 
-0.00284495, 0, 0.00331832, 0.00283554, 0.00141476, -0.00141476, 
-0.00188947, 0.00141743, -0.00236351, 0.00236351, 0.00235794, 
0.00235239, -0.000940292, -0.0014121, -0.00283019, 0.000472255, 
0.000472032, 0.000471809, -0.0014161, 0.0014161, -0.000943842, 
0.000472032, -0.000944287, -0.00094518, -0.00189304, -0.000473821, 
-0.000474046, 0.00331361, -0.000472701, -0.000946074, 0.00141878, 
-0.000945627, -0.00189394, -0.00189753, -0.0057143, -0.00143369, 
-0.00383326, 0.00143919, 0.000479272, -0.00191847, -0.000480192, 
0.000960154, 0.000479731, 0, 0.000479501, 0.000958313, -0.00383878, 
-0.00240674, 0.000963391, 0.000962464, -0.00192586, 0.000481812, 
-0.00241138, -0.00144963)

for(i in 6:0)
	cat(i, try(fitdist(x2*10^i, "cauchy", method="mle")$estimate, silent=TRUE), "\n")


# (15) Fit of a lognormal distribution on acute toxicity values of endosulfan for
# nonarthropod invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5 percent quantile value of 
# the fitted distribution, what is called the 5 percent hazardous concentration (HC5)
# in ecotoxicology, with its two-sided 95 percent confidence interval calculated by 
# parametric bootstrap
#
data(endosulfan)
ATV <- subset(endosulfan, group == "NonArthroInvert")$ATV
log10ATV <- log10(subset(endosulfan, group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV, "norm")

quantile(fln, probs = 0.05, bootstrap=TRUE, 
	bootstrap.arg = list(bootmethod = "param", niter = 501))

quantile(fln, probs = c(0.05,0.1,0.2), bootstrap = TRUE, CI.type = "greater",
	bootstrap.arg = list(bootmethod = "param", niter = 501))

quantile(fln, probs = 0.05, bootstrap=TRUE, 
	bootstrap.arg = list(bootmethod = "nonparam", niter = 501))



# (16) uniform distribution
#

x3 <- runif(1000, 1/2, 5/2)
fitdist(x3, "unif", method="mle")
# mledist(x3, "unif")
fitdist(x3, "unif", method="mge", gof="CvM")
# mgedist(x3, "unif")
fitdist(x3, "unif", method="qme", probs=c(1/10, 9/10))




# (17) uniform distribution
#
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

data(danishuni)
fitdist(danishuni$Loss, "gumbel", start=list(a=5, b=10))


# (18) check the 'start' argument
#

#create the sample
mysample <- rexp(100, 5)
mystart2 <- list(rate2=8)
mystart3 <- list(8)

try( fitdist(mysample, dexp, start= mystart2, method="mle") ) 
try( fitdist(mysample, dexp, start= mystart3, method="mle") ) 

try( fitdist(mysample, dexp, start= mystart2, method="mme") ) 
try( fitdist(mysample, dexp, start= mystart3, method="mme") ) 

try( fitdist(mysample, dexp, start= mystart2, method="qme", probs=1/2) ) 
try( fitdist(mysample, dexp, start= mystart3, method="qme", probs=1/2) ) 

try( fitdist(mysample, dexp, start= mystart2, method="mge", gof="AD") ) 
try( fitdist(mysample, dexp, start= mystart3, method="mge", gof="AD") ) 

