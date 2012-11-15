library(fitdistrplus)



# (1) basic fit of a normal distribution with maximum likelihood estimation
#

x1<-c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
mledist(x1,"norm")

# (2) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view dedicated to probability distributions

dgumbel<-function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
mledist(x1,"gumbel",start=list(a=10,b=5))
mledist(x1,"gumbel",start=list(a=10), fix.arg=list(b=3))

# (3) fit a discrete distribution (Poisson)
#

x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
mledist(x2,"pois")
mledist(x2,"nbinom")

# (4) fit a finite-support distribution (beta)
#

x3<-c(0.80,0.72,0.88,0.84,0.38,0.64,0.69,0.48,0.73,0.58,0.81,
0.83,0.71,0.75,0.59)
mledist(x3,"beta")


# (5) fit frequency distributions on USArrests dataset.
#

x4 <- USArrests$Assault
mledist(x4, "pois")
mledist(x4, "nbinom")

# (6) fit a continuous distribution (Gumbel) to censored data.
#

d1<-data.frame(
left=c(1.73,1.51,0.77,1.96,1.96,-1.4,-1.4,NA,-0.11,0.55,0.41,
2.56,NA,-0.53,0.63,-1.4,-1.4,-1.4,NA,0.13),
right=c(1.73,1.51,0.77,1.96,1.96,0,-0.7,-1.4,-0.11,0.55,0.41,
2.56,-1.4,-0.53,0.63,0,-0.7,NA,-1.4,0.13))
mledist(d1,"norm")

dgumbel<-function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel<-function(q,a,b) exp(-exp((a-q)/b))
mledist(d1,"gumbel",start=list(a=0,b=2),optim.method="Nelder-Mead")


# (7) scaling problem
#

x <- c(-0.00707717, -0.000947418, -0.00189753, 
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


#only i == 0, no scaling, should not converge.
for(i in 6:0)
cat(i, try(mledist(x*10^i, "cauchy")$estimate, silent=TRUE), "\n")


# (8) normal mixture
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
#MLE fit
fit1 <- mledist(x, "norm2", start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2), 
	lower=c(0, 0, 0, 0, 0))






# (9) fit a Pareto distribution
#

if(any(installed.packages()[,"Package"] == "actuar"))
{
    require(actuar)
#simulate a sample
    x4 <- rpareto(1000, 6, 2)
	
#fit
    mledist(x4, "pareto", start=c(10, 10), lower=1, upper=Inf)
	
}



# (10) custom optim for exponential distribution
#
if(any(installed.packages()[,"Package"] == "rgenoud"))
{
	

mysample <- rexp(1000, 5)
mystart <- 8

fNM <- mledist(mysample, "exp", optim.method="Nelder-Mead") 
fBFGS <- mledist(mysample, "exp", optim.method="BFGS") 
fLBFGSB <- mledist(mysample, "exp", optim.method="L-BFGS-B", lower=0) 
fSANN <- mledist(mysample, "exp", optim.method="SANN") 
fCG <- mledist(mysample, "exp", optim.method="CG") 

#the warning tell us to use optimise...

#to meet the standard 'fn' argument and specific name arguments, we wrap optimize,
myoptimize <- function(fn, par, ...)
{
	res <- optimize(f=fn, ..., maximum=FALSE)	
	c(res, convergence=0, value=res$objective, par=res$minimum, hessian=NA)
}

foptimize <- mledist(mysample, "exp", start=mystart, custom.optim=myoptimize, interval=c(0, 100))


library(rgenoud)

#wrap genoud function rgenoud package
mygenoud <- function(fn, par, ...) 
{
	res <- genoud(fn, starting.values=par, ...)        
	c(res, convergence=0)       
}


fgenoud <- mledist(mysample, "exp", start=mystart, custom.optim= mygenoud, nvars=1,    
        Domains=cbind(0, 10), boundary.enforcement=1, 
        hessian=TRUE, print.level=0)


c(NM=fNM$estimate,
BFGS=fBFGS$estimate,
LBFGSB=fLBFGSB$estimate,
SANN=fSANN$estimate,
CG=fCG$estimate,
optimize=foptimize$estimate,
fgenoud=fgenoud$estimate)

}


# (10) custom optim for gamma distribution
#
if(any(installed.packages()[,"Package"] == "rgenoud"))
{
	

mysample <- rgamma(1000, 5, 3)
mystart <- c(10, 10)

fNM <- mledist(mysample, "gamma", optim.method="Nelder-Mead") 
fBFGS <- mledist(mysample, "gamma", optim.method="BFGS") 
fLBFGSB <- mledist(mysample, "gamma", optim.method="L-BFGS-B", lower=0) 
fSANN <- mledist(mysample, "gamma", optim.method="SANN") 
fCG <- mledist(mysample, "gamma", optim.method="CG", control=list(maxit=1000)) 

fgenoud <- mledist(mysample, "gamma", start=mystart, custom.optim= mygenoud, nvars=2,    
        Domains=cbind(c(0,0), c(100,100)), boundary.enforcement=1, 
        hessian=TRUE, print.level=0)

cbind(NM=fNM$estimate,
BFGS=fBFGS$estimate,
LBFGSB=fLBFGSB$estimate,
SANN=fSANN$estimate,
CG=fCG$estimate,
fgenoud=fgenoud$estimate)


data(groundbeef)

fNM <- mledist(groundbeef$serving, "gamma", optim.method="Nelder-Mead") 
fBFGS <- mledist(groundbeef$serving, "gamma", optim.method="BFGS") 
fLBFGSB <- mledist(groundbeef$serving, "gamma", optim.method="L-BFGS-B", lower=0) 
fSANN <- mledist(groundbeef$serving, "gamma", optim.method="SANN") 
fCG <- mledist(groundbeef$serving, "gamma", optim.method="CG", control=list(maxit=10000)) 

fgenoud <- mledist(groundbeef$serving, "gamma", start=mystart, 
		custom.optim= mygenoud, nvars=2,    
        Domains=cbind(c(0,0), c(100,100)), boundary.enforcement=1, 
        hessian=TRUE, print.level=0)

cbind(NM=fNM$estimate,
BFGS=fBFGS$estimate,
LBFGSB=fLBFGSB$estimate,
SANN=fSANN$estimate,
CG=fCG$estimate,
fgenoud=fgenoud$estimate)


}


