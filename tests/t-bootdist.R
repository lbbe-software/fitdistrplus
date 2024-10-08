require("fitdistrplus")
#We choose a low number of bootstrap replicates in order to satisfy CRAN running times constraint.
#For practical application, we recommend to use nbboot=501 or nbboot=1001.

nbboot <- 1001
nbboot <- 11

nsample <- 100
nsample <- 10

visualize <- FALSE # TRUE for manual tests with visualization of results

# (1) Fit of a gamma distribution to serving size data
# using default method (maximum likelihood estimation)
# followed by parametric bootstrap
#
data(groundbeef)
serving <- groundbeef$serving
f1 <- fitdist(serving, "gamma")
b1 <- bootdist(f1, niter=nbboot, silent=TRUE)
b1 <- bootdist(f1, niter=nbboot, silent=FALSE)

print(lapply(b1, head))
plot(b1)
summary(b1)

# (2) new plot arguments
#for new graph functions
f1 <- fitdist(rgamma(nsample, 2, 3), "gamma")
b1 <- bootdist(f1, niter=nbboot, silent=TRUE)
plot(b1)
plot(b1, trueval = c(2, 3))
plot(b1, enhance=TRUE)
plot(b1, enhance=TRUE, trueval = c(2, 3))
plot(b1, enhance=TRUE, rampcol=c("blue", "green"), nbgrid=15, nbcol=15)

if(any(installed.packages()[, "Package"] == "actuar") && visualize)
{
  require("actuar")
  set.seed(123)
  f1 <- fitdist(rburr(nsample, 2, 3, 1), "burr", start=list(shape1=10, shape2=10, rate=1))
  b1 <- bootdist(f1, niter=nbboot, silent=TRUE)
  plot(b1)
  plot(b1, trueval = c(2, 3, 1))
  plot(b1, enhance=TRUE)
  plot(b1, enhance=TRUE, trueval = c(2, 3, 1))
}


# (3) estimation of the rate of a gamma distribution 
# by maximum likelihood with the shape fixed at 4 using the argument fix.arg
# followed by parametric bootstrap
#
f1c  <- fitdist(serving, "gamma", start=list(rate=0.1), fix.arg=list(shape=4))
b1c  <-  bootdist(f1c, niter=nbboot)
summary(b1c)

# (4) fit of a gamma distribution to serving size data 
# by quantile matching estimation (in this example matching 
# first and third quartiles) followed by parametric bootstrap
#
f1d  <-  fitdist(serving, "gamma", method="qme", probs=c(0.25, 0.75))
b1d  <-  bootdist(f1d, niter=nbboot)
summary(b1d)

# (5) fit of a gamma distribution with control of the optimization
# method,  followed by parametric bootstrap
#
if(visualize) { # check ERROR on aarch64-apple-darwin20.4.0 (64-bit) (2021/05/12)
  set.seed(1234)
  f1e <- fitdist(serving, "gamma", optim.method="L-BFGS-B", lower=c(0, 0))
  b1e <- bootdist(f1e, niter=nbboot)
  summary(b1e)
}

# (6) fit of a discrete distribution by matching moment estimation 
# (using a closed formula) followed by parametric bootstrap
#
set.seed(1234)
x2 <- rpois(nsample, lambda = 5)
f2 <- fitdist(x2, "pois", method="mme")
b2 <- bootdist(f2, niter=nbboot)
plot(b2,pch=16)
summary(b2)

# (7) Fit of a uniform distribution using the Cramer-von Mises distance
# followed by parametric bootstrap 
# 
if(visualize)
{
  x3  <-  runif(nsample, min=5, max=10)
  f3  <-  fitdist(x3, "unif", method="mge", gof="CvM")
  b3  <-  bootdist(f3, bootmethod="param", niter=nbboot)
  summary(b3)
  plot(b3)
}

# (9) fit of a Weibull distribution to serving size data by maximum likelihood 
# estimation or by quantile matching estimation (in this example matching 
# first and third quartiles) followed by parametric bootstrap
#
fWmle <- fitdist(serving, "weibull")
bWmle <- bootdist(fWmle, niter=nbboot)
summary(bWmle)
quantile(bWmle, probs=c(0.25, 0.75))

fWqme <- fitdist(serving, "weibull", method="qme", probs=c(0.25, 0.75))
bWqme <- bootdist(fWqme, niter=nbboot)
summary(bWqme)
quantile(bWqme, probs=c(0.25, 0.75))


# (10) Fit of a Pareto distribution by numerical moment matching estimation
# followed by parametric bootstrap
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LONG TO RUN !!!!!!!!!!!!!!!!!!!!!!!!
#
if (visualize) 
{
  if(any(installed.packages()[, "Package"] == "actuar"))
  {
    require("actuar")
    #simulate a sample
    x4 <- rpareto(nsample,  6,  2)
    memp <- function(x,  order)
      ifelse(order == 1,  mean(x),  sum(x^order)/length(x))
    
    f4 <- fitdist(x4,  "pareto",  "mme",  order=1:2,  
                  start=list(shape=10,  scale=10),  
                  lower=1,  memp="memp",  upper=50)
    
    b4 <- bootdist(f4,  niter=nbboot)
    summary(b4)
    
    b4npar <- bootdist(f4,  niter=nbboot,  bootmethod="nonparam")
    summary(b4npar)
  }
}

# (11) Fit of a Burr distribution (3 parameters) using MLE 
# followed by parametric boostrap
# !!!!!!!!!!!!!!!! LONG TO RUN !!!!!!!!!!!!!!!!!!
#
if (visualize)
{
  if(any(installed.packages()[, "Package"] == "actuar"))
  {
    require("actuar")
    data(danishuni)
    
    fdan <- fitdist(danishuni$Loss, "burr", method="mle", 
                    start=list(shape1=5, shape2=5, rate=10), lower=0+1e-1, control=list(trace=0))
    bdan <- bootdist(fdan,  bootmethod="param", niter=nbboot)
    summary(bdan)
    plot(bdan)
    cdfcomp(fdan, xlogscale=TRUE)
  }
}

# (12) Fit of a Triangular distribution (3 parameters) using MLE 
# followed by parametric boostrap, with crashes of optim
# 
if(any(installed.packages()[, "Package"] == "mc2d"))
{
  require("mc2d")
  set.seed(1234)
  x4 <- rtriang(100,min=0,mode=4,max=20) # nsample not used : does not converge if the sample is too small
  fit4t<-fitdist(x4,dtriang,start=list(min=0,mode=4,max=20))
  summary(fit4t)
  b4t<-bootdist(fit4t,niter=nbboot) 
  b4t
  plot(b4t)
  summary(b4t)
  quantile(b4t)
  
}

# (13) Fit of a Pareto and a Burr distribution, with bootstrap on the Burr distribution
#
#
if(visualize)
{
  data(endosulfan)
  ATV <-endosulfan$ATV
  plotdist(ATV)
  descdist(ATV,boot=nbboot)
  
  fln <- fitdist(ATV, "lnorm")
  summary(fln)
  gofstat(fln)
  
  # use of plotdist to find good reasonable initial values for parameters
  plotdist(ATV, "pareto", para=list(shape=1, scale=500))
  fP <- fitdist(ATV, "pareto", start=list(shape=1, scale=500))
  summary(fP)
  gofstat(fP)
  
  # definition of the initial values from the fit of the Pareto
  # as the Burr distribution is the Pareto when shape2 == 1
  fB <- fitdist(ATV, "burr", start=list(shape1=0.3, shape2=1, rate=1))
  summary(fB)
  gofstat(fB)
  
  cdfcomp(list(fln,fP,fB),xlogscale=TRUE)
  qqcomp(list(fln,fP,fB),xlogscale=TRUE,ylogscale=TRUE)
  ppcomp(list(fln,fP,fB),xlogscale=TRUE,ylogscale=TRUE)
  denscomp(list(fln,fP,fB)) # without great interest as hist does accept argument log="x"
  
  # comparison of HC5 values (5 percent quantiles)
  quantile(fln,probs=0.05)
  quantile(fP,probs=0.05)
  quantile(fB,probs=0.05)
  
  # bootstrap for the Burr distribution
  bfB <- bootdist(fB,niter=nbboot)
  plot(bfB)
}

# (14) relevant example for zero modified geometric distribution
#
dzmgeom <- function(x, p1, p2)
{
  p1 * (x == 0) + (1-p1)*dgeom(x-1, p2)
}
pzmgeom <- function(q, p1, p2)
{
  p1 * (q >= 0) + (1-p1)*pgeom(q-1, p2)
}
rzmgeom <- function(n, p1, p2)
{
  u <- rbinom(n, 1, 1-p1) #prob to get zero is p1
  u[u != 0] <- rgeom(sum(u != 0), p2)+1
  u
}

x2 <- rzmgeom(nsample, 1/2, 1/10)

f2 <- fitdist(x2, "zmgeom", method="mle", fix.arg=function(x) list(p1=mean(x == 0)), start=list(p2=1/2))
b2 <- bootdist(f2, niter=nbboot)
plot(b2)

f3 <- fitdist(x2, "zmgeom", method="mle", start=list(p1=1/2, p2=1/2))
b3 <- bootdist(f3, niter=nbboot)
plot(b3, enhance=TRUE)

# (15) does fixing p1 reduce bias of estimating p2?
summary(b2$estim[, "p2"] - 1/10)
summary(b3$estim[, "p2"] - 1/10)

par(mfrow=c(1, 2))
hist(b2$estim[, "p2"] - 1/10, breaks=100, xlim=c(-.015, .015))
hist(b3$estim[, "p2"] - 1/10, breaks=100, xlim=c(-.015, .015))
par(mfrow=c(1, 1))

# (16) efficiency of parallel operation
if (visualize)
{
  niter <- 1001 
  data(groundbeef)
  serving <- groundbeef$serving
  f1 <- fitdist(serving, "gamma")
  
  alltime <- matrix(NA, 9, 5)
  colnames(alltime) <- c("user.self",  "sys.self",   "elapsed",    "user.child", "sys.child" )
  rownames(alltime) <- c("base R", paste("snow", 1:4), paste("multicore", 1:4))
  
  alltime[1,] <- system.time(res <- bootdist(f1, niter = niter))
  
  for (cli in 1:4)
  {
    cat("\nnb cluster", cli, "\n")
    #ptm <- proc.time()
    alltime[cli+1,] <- system.time(res <- bootdist(f1, niter = niter, parallel = "snow", ncpus = cli))
    print(summary(res))
    #print(proc.time() - ptm)
    
  }
  # not available on Windows
  if(.Platform$OS.type == "unix") 
    for (cli in 1:4)
    {
      cat("\nnb cluster", cli, "\n")
      #ptm <- proc.time()
      alltime[cli+5,] <- system.time(res <- bootdist(f1, niter = niter, parallel = "multicore", ncpus = cli))
      print(summary(res))
      #print(proc.time() - ptm)
      
    }
  
  alltime
}

# (17) bootdist with weights (not yet available, test of error message)
#
x <- rpois(nsample, 10)
xtab <- table(x)
xval <- sort(unique(x))
(f1 <- fitdist(x, "pois"))
(f2 <- fitdist(xval, "pois", weights = xtab))
summary(bootdist(f1, niter = nbboot))
try(summary(bootdist(f2, niter = nbboot))) # not yet developed

# (18) density of bootdist()
#
x <- rlnorm(1e3)
b0 <- bootdist(fitdist(x, "lnorm"), niter = 50)
b1 <- bootdist(fitdist(x, "lnorm"), niter = 100)
b2 <- bootdist(fitdist(x, "lnorm"), niter = 200)

#d1 <- fitdistrplus:::density.bootdist(b0, b1, b2)
d1 <- density(b0, b1, b2)
str(d1)
plot(d1)
print(d1)

d1 <- density(b1)
str(d1)
plot(d1)
print(d1)