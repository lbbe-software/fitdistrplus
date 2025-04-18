require("fitdistrplus")

data("groundbeef")
serving <- groundbeef$serving 

# (1) basic fit of a normal distribution with maximum likelihood estimation
#

set.seed(1234)
x1 <- rnorm(n=100)
mledist(x1,"norm")

#fitted coef around -0.1567617  0.9993707, fitted loglik -141.8309  

# (2) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view dedicated to probability distributions

dgumbel <- function(x,a,b)
{
  cat(a, b, "\n")
  res <- 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
  if(any(is.infinite(log(res))))
    stop(log(res))
  res
}
mledist(x1,"gumbel",start=list(a=10,b=5), silent=TRUE)
mledist(x1,"gumbel",start=list(a=10,b=5), silent=FALSE, control=list(trace=1), lower=0)

#fitted coef around -0.6267919  0.8564855, fitted loglik -139.3812


# (3) fit a discrete distribution (Poisson)
#

set.seed(1234)
x2 <- rpois(n=30,lambda = 2)
mledist(x2,"pois")

#fitted coef around 1.7, fitted loglik -46.18434

# (4) fit a finite-support distribution (beta)
#

set.seed(1234)
x3 <- rbeta(n=100,shape1=5, shape2=10)
mledist(x3,"beta")

#fitted coef around 4.859798 10.918841, fitted loglik 78.33052

# (5) fit frequency distributions on USArrests dataset.
#

x4 <- USArrests$Assault
mledist(x4, "pois", silent=TRUE)
mledist(x4, "pois", silent=FALSE)
mledist(x4, "nbinom")

#fitted coef around 3.822579 170.747853, fitted loglik -290.3297


# (6) scaling problem
# the simulated dataset (below) has particularly small values, hence without scaling (10^0),
# the optimization raises an error. The for loop shows how scaling by 10^i
# for i=1,...,6 makes the fitting procedure work correctly.

set.seed(1234)
x2 <- rnorm(100, 1e-4, 2e-4)
for(i in 6:0)
  cat(i, try(mledist(x*10^i, "cauchy")$estimate, silent=TRUE), "\n")


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
dnorm2 <- function(x, w, m1, s1, m2, s2)
  w*dnorm(x, m1, s1) + (1-w)*dnorm(x, m2, s2)
#distribution function		
pnorm2 <- function(q, w, m1, s1, m2, s2)
  w*pnorm(q, m1, s1) + (1-w)*pnorm(q, m2, s2)		
#numerically-approximated quantile function
qnorm2 <- function(p, w, m1, s1, m2, s2)
{
  L2 <- function(x, prob)
    (prob - pnorm2(x, w, m1, s1, m2, s2))^2	
  sapply(p, function(pr) optimize(L2, c(-20, 30), prob=pr)$minimum)
}	


#basic normal distribution
x <- c(rnorm(1000, 5),  rnorm(1000, 10))
#MLE fit
fit1 <- mledist(x, "norm2", start=list(w=1/3, m1=4, s1=2, m2=8, s2=2), 
                lower=c(0, 0, 0, 0, 0))
fit1

#fitted coef around 0.5003298  4.9842719  0.9909527 10.0296973  1.0024444 , fitted loglik -4185.114

# (9) fit a Pareto and log-logistic distributions
#

if(any(installed.packages()[,"Package"] == "actuar"))
{
  require(actuar)
  
  #simulate a sample
  x4 <- rpareto(1000, 6, 2)
  
  #fit
  mledist(x4, "pareto", start=list(shape=10, scale=10), lower=1, upper=Inf)
  
  #fitted coef around 6.334596 2.110411, fitted loglik -58.73242
  
  #simulate a sample
  x4 <- rllogis(1000, 6, 2)
  
  #fit
  mledist(x4, "llogis", start=list(shape=10, rate=2), lower=1, upper=Inf)
  
  #fitted coef around 6.150001 2.005405, fitted loglik 511.873
}



# (10) custom optim for exponential distribution
#
if(any(installed.packages()[,"Package"] == "rgenoud") && FALSE)
{
  
  
  mysample <- rexp(1000, 5)
  mystart <- list(rate=8)
  
  fNM <- mledist(mysample, "exp", optim.method="Nelder-Mead") 
  fBFGS <- mledist(mysample, "exp", optim.method="BFGS") 
  fLBFGSB <- mledist(mysample, "exp", optim.method="L-BFGS-B", lower=0) 
  fSANN <- mledist(mysample, "exp", optim.method="SANN") 
  fCG <- try(mledist(mysample, "exp", optim.method="CG") )
  if(inherits(fCG, "try-error"))
    fCG <- list(estimate=NA)
  
  #the warning tell us to use optimise...
  
  #to meet the standard 'fn' argument and specific name arguments, we wrap optimize,
  myoptimize <- function(fn, par, ...)
  {
    res <- optimize(f=fn, ..., maximum=FALSE)	
    c(res, convergence=0, value=res$objective, par=res$minimum, hessian=NA)
  }
  
  foptimize <- mledist(mysample, "exp", start=mystart, custom.optim=myoptimize, interval=c(0, 100))
  
  
  require("rgenoud")
  
  #wrap genoud function rgenoud package
  mygenoud <- function(fn, par, ...) 
  {
    res <- genoud(fn, starting.values=par, ...)        
    c(res, convergence=0, counts=NULL)       
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


# (11) custom optim for gamma distribution
#
if(any(installed.packages()[,"Package"] == "rgenoud") && FALSE)
{
  
  
  mysample <- rgamma(1000, 5, 3)
  mystart <- c(shape=10, rate=10)
  
  fNM <- mledist(mysample, "gamma", optim.method="Nelder-Mead") 
  fBFGS <- mledist(mysample, "gamma", optim.method="BFGS") 
  fLBFGSB <- mledist(mysample, "gamma", optim.method="L-BFGS-B", lower=0) 
  fSANN <- mledist(mysample, "gamma", optim.method="SANN") 
  fCG <- try( mledist(mysample, "gamma", optim.method="CG", control=list(maxit=1000)) )
  if(inherits(fCG, "try-error"))
    fCG <- list(estimate=NA)
  
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
  fCG <- try( mledist(groundbeef$serving, "gamma", optim.method="CG", control=list(maxit=10000)) )
  if(inherits(fCG, "try-error"))
    fCG <- list(estimate=NA)
  
  fgenoud <- mledist(groundbeef$serving, "gamma", start=list(shape=4, rate=1),
                     custom.optim= mygenoud, nvars=2, max.generations=10,
                     Domains=cbind(c(0,0), c(10,10)), boundary.enforcement=1, 
                     hessian=TRUE, print.level=0, P9=10)
  
  cbind(NM=fNM$estimate,
        BFGS=fBFGS$estimate,
        LBFGSB=fLBFGSB$estimate,
        SANN=fSANN$estimate,
        CG=fCG$estimate,
        fgenoud=fgenoud$estimate)
  
  
}



# (12) test error messages
#

dnorm2 <- function(x, a)
  "NA"
x <- rexp(10)

#should get a one-line error 
res <- mledist(x, "norm2", start=list(a=1))
#as in 
attr(try(log("a"), silent=TRUE), "condition")


try(mledist(c(serving, "a"), "gamma"))
try(mledist(c(serving, NA), "gamma"))
try(mledist(c(serving, Inf), "gamma"))
try(mledist(c(serving, -Inf), "gamma"))
try(mledist(c(serving, NaN), "gamma"))



# (13) weighted MLE
#
n <- 1e6
n <- 1e2
x <- rpois(n, 10)
xtab <- table(x)
xval <- sort(unique(x))
f1 <- mledist(x, "pois", start=list(lambda=mean(x)), optim.method="Brent", lower=0, 
              upper=100, control=list(trace=1))
f2 <- mledist(xval, "pois", weights=xtab, start=list(lambda=mean(x)))

f1$estimate
f2$estimate #should be identical

#fitted coef around 9.74, fitted loglik -246.8505

#test discrete distrib
f2 <- try(mledist(xval, "pois", weights=1:length(xval), start=list(lambda=mean(x))))
#test non integer weights
f2 <- try(mledist(xval, "pois", weights=rep(1/3, length(xval)), start=list(lambda=mean(x))))
f2 <- try(mledist(1:10, "pois", weights=c(rep(1, 9), 1.001), start=list(lambda=mean(x))))
f2 <- try(mledist(1:10, "pois", weights=c(rep(1, 9), 1.0000001), start=list(lambda=mean(x))))

# (14) no convergence
#
n <- 1e2
x <- c(rep(0, n), rpois(n, 10), rpois(n, 50))

mledist(x, "pois", optim.method="Nelder-Mead", control=list(maxit=10))


# (15) basic fit of a normal distribution with new fix.arg/start.arg
#

set.seed(1234)
x1 <- rnorm(n=100)

#correct usage
mledist(x1,"norm")
mledist(x1,"norm", start=function(x) list(mean=0, sd=1))
mledist(x1,"norm", fix.arg=function(x) list(mean=mean(x)))
mledist(x1,"norm", fix.arg=list(mean=1/2))
mledist(x1,"norm", fix.arg=list(mean=1/2), start=list(sd=1))
mledist(x1,"norm", fix.arg=function(x) list(mean=0), start=list(sd=1))
mledist(x1,"norm", fix.arg=function(x) list(mean=mean(x)), start=list(sd=1))

#wrong usage (see text message in util-checkparam.R)
try( mledist(x1,"norm", start=list(a=1/2)) ) #t3
try( mledist(x1,"norm", start=function(x) list(a=0, b=1)) ) #t3
try( mledist(x1,"norm", fix.arg=list(a=1/2)) ) #t4
try( mledist(x1,"norm", fix.arg=function(x) list(a=0), start=list(sd=1)) ) #t4
try( mledist(x1,"norm", start=matrix(1/2)) ) #t1
try( mledist(x1,"norm", fix.arg=matrix(1/2)) ) #t0
try( mledist(x1,"norm", fix.arg=matrix(1/2), start=matrix(1/2)) ) #t2
try( mledist(x1,"norm", fix.arg=function(x) list(mean=mean(x), sd=2), start=list(sd=1)) ) #t5
dabcnorm <- function(x, mean, sd) 1
try( mledist(x1,"abcnorm", fix.arg=function(x) list(mean=mean(x))) ) #t8



# (16) relevant example for zero modified geometric distribution
#
dzmgeom <- function(x, p1, p2)
{
  p1 * (x == 0) + (1-p1)*dgeom(x-1, p2)
}
rzmgeom <- function(n, p1, p2)
{
  u <- rbinom(n, 1, 1-p1) #prob to get zero is p1
  u[u != 0] <- rgeom(sum(u != 0), p2)+1
  u
}
# check
# dzmgeom(0:5, 1/2, 1/10)

x2 <- rzmgeom(100, 1/2, 1/10)
x3 <- rzmgeom(100, 1/3, 1/10)
x4 <- rzmgeom(100, 1/4, 1/10)
table(x2)
#this is the MLE which converges almost surely and in distribution to the true value.
initp1 <- function(x) list(p1=mean(x == 0))

mledist(x2, "zmgeom", fix.arg=initp1, start=list(p2=1/2))[c("estimate", "fix.arg")]
mledist(x3, "zmgeom", fix.arg=initp1, start=list(p2=1/2))[c("estimate", "fix.arg")]
mledist(x4, "zmgeom", fix.arg=initp1, start=list(p2=1/2))[c("estimate", "fix.arg")]


#fitted coef around 0.09044335, fitted loglik -298.5162


# (17) test the component optim.message
x <- rnorm(100)
#change parameter to obtain unsuccessful convergence
mledist(x, "norm", control=list(maxit=2), start=list(mean=1e5, sd=1), optim.method="L-BFGS-B", lower=0)


# (18) management of bounds in optim/constrOptim
x <- rexp(100)
mledist(x, "exp") #optim, BFGS
mledist(x, "exp", optim.method="Brent", lower=0, upper=100) #optim, Brent
mledist(x, "exp", optim.method="Nelder-Mead") #optim, Nelder-Mead
mledist(x, "exp", lower=0, optim.method="Nelder-Mead") #constrOptim, Nelder-Mead
mledist(x, "exp", lower=0, optim.method="BFGS") #optim, L-BFGS-B


x <- rbeta(100, 3/2, 7/3)
mledist(x, "beta", optim.method="Nelder") #optim, Nelder-Mead
mledist(x, "beta", lower=0, optim.method="Nelder-Mead") #constrOptim, Nelder-Mead
#as the result of optim(c(-1.2,1), fr, method = "Nelder-Mead", hessian=TRUE, gr=NULL, lower=-Inf, upper=Inf) from optim() example

mledist(x, "beta", lower=0, optim.method="BFGS") #optim, L-BFGS-B


# (19) large sample size issue

if(FALSE)
{
  set.seed(123)
  obs <- rlnorm(1e7, 3, 2)
  for(i in 2:7)
    cat(i, try(mledist(obs[1:10^i], "lnorm", control=list(trace=0, REPORT=1))$estimate, silent=TRUE), "\n")
  
  # 2 3.143734 1.819549 
  # 3 3.023688 1.955416 
  # 4 2.995646 1.993392 
  # 5 3.002682 2.000394 
  # 6 2.999036 1.999887 
  # 7 3.000222 1.999981 
}

# (20) threshold parameter leading to infinite log-density

set.seed(123)
dsexp <- function(x, rate, shift)
  dexp(x-shift, rate=rate)
psexp <- function(x, rate, shift)
  pexp(x-shift, rate=rate)
rsexp <- function(n, rate, shift)
  rexp(n, rate=rate)+shift
x <- rsexp(1000, 1/4, 1)
sum(log(dsexp(x, 1, 1.995)))
mledist(x, "sexp", start=list(rate=1, shift=0), upper= c(Inf, min(x)), control=list(trace=1, REPORT=1), optim="L-BFGS-B")
mledist(x, "sexp", start=list(rate=1, shift=0), upper= c(Inf, min(x)), control=list(trace=1, REPORT=1))

#fitted coef is 0.243 1.003 , fitted loglik -2415, fitted hessian is        
#       rate     shift
# rate   1000 -1.00e+03
# shift -1000  2.14e-05




