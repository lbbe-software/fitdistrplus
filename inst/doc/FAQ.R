## ----setup, echo=FALSE, message=FALSE, warning=FALSE---------------------
require(fitdistrplus)
set.seed(1234)
options(digits = 3)

## ---- eval=FALSE---------------------------------------------------------
#  dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
#  pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
#  qgumbel <- function(p, a, b) a-b*log(-log(p))
#  data(groundbeef)
#  fitgumbel <- fitdist(groundbeef$serving, "gumbel", start=list(a=10, b=10))

## ---- eval=FALSE---------------------------------------------------------
#  dzmgeom <- function(x, p1, p2) p1 * (x == 0) + (1-p1)*dgeom(x-1, p2)
#  pzmgeom <- function(q, p1, p2) p1 * (q >= 0) + (1-p1)*pgeom(q-1, p2)
#  rzmgeom <- function(n, p1, p2)
#  {
#    u <- rbinom(n, 1, 1-p1) #prob to get zero is p1
#    u[u != 0] <- rgeom(sum(u != 0), p2)+1
#    u
#  }
#  x2 <- rzmgeom(1000, 1/2, 1/10)
#  fitdist(x2, "zmgeom", start=list(p1=1/2, p2=1/2))

## ---- message=FALSE------------------------------------------------------
data("endosulfan")
library("actuar")
fendo.B <- fitdist(endosulfan$ATV, "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1))
summary(fendo.B)

## ---- fig.height=3, fig.width=6------------------------------------------
x3 <- rlnorm(1000)
f1 <- fitdist(x3, "lnorm", method="mle") 
f2 <- fitdist(x3, "lnorm", method="mme")
par(mfrow=1:2)
cdfcomp(list(f1, f2), do.points=FALSE, xlogscale = TRUE, main = "CDF plot")
denscomp(list(f1, f2), demp=TRUE, main = "Density plot")

## ------------------------------------------------------------------------
c("E(X) by MME"=as.numeric(exp(f2$estimate["meanlog"]+f2$estimate["sdlog"]^2/2)), 
	"E(X) by MLE"=as.numeric(exp(f1$estimate["meanlog"]+f1$estimate["sdlog"]^2/2)), 
	"empirical"=mean(x3))
c("Var(X) by MME"=as.numeric(exp(2*f2$estimate["meanlog"]+f2$estimate["sdlog"]^2)*(exp(f2$estimate["sdlog"]^2)-1)), 
	"Var(X) by MLE"=as.numeric(exp(2*f1$estimate["meanlog"]+f1$estimate["sdlog"]^2)*(exp(f1$estimate["sdlog"]^2)-1)), 
	"empirical"=var(x3))

## ------------------------------------------------------------------------
set.seed(1234)
x <- rnorm(100, mean = 1, sd = 0.5)
(try(fitdist(x, "exp")))

## ------------------------------------------------------------------------
fitdist(x[x >= 0], "exp")
fitdist(x - min(x), "exp")

## ------------------------------------------------------------------------
set.seed(1234)
x <- rnorm(100, mean = 0.5, sd = 0.25)
(try(fitdist(x, "beta")))

## ------------------------------------------------------------------------
fitdist(x[x > 0 & x < 1], "beta")
fitdist((x - min(x)*1.01) / (max(x) * 1.01 - min(x) * 1.01), "beta")

## ---- fig.height=3, fig.width=6------------------------------------------
set.seed(1234)
x <- rgamma(n = 100, shape = 2, scale = 1)
# fit of the good distribution
fgamma <- fitdist(x, "gamma")
# fit of a bad distribution
fexp <- fitdist(x, "exp")
g <- gofstat(list(fgamma, fexp), fitnames = c("gamma", "exp"))
denscomp(list(fgamma, fexp), legendtext = c("gamma", "exp"))
# results of the tests
## chi square test (with corresponding table with theoretical and observed counts)
g$chisqpvalue
g$chisqtable
## Anderson-Darling test
g$adtest
## Cramer von  Mises test
g$cvmtest
## Kolmogorov-Smirnov test
g$kstest

## ---- fig.height=3, fig.width=6------------------------------------------
set.seed(1234)
x1 <- rpois(n = 100, lambda = 100)
f1 <- fitdist(x1, "norm")
g1 <- gofstat(f1)
g1$kstest

x2 <- rpois(n = 10000, lambda = 100)
f2 <- fitdist(x2, "norm")
g2 <- gofstat(f2)
g2$kstest

par(mfrow=1:2)
denscomp(f1, demp = TRUE, addlegend = FALSE, main = "small sample")
denscomp(f2, demp = TRUE, addlegend = FALSE, main = "big sample")


## ---- fig.height=3, fig.width=6------------------------------------------
set.seed(1234)
x3 <- rpois(n = 500, lambda = 1)
f3 <- fitdist(x3, "norm")
g3 <- gofstat(f3)
g3$kstest

x4 <- rpois(n = 50, lambda = 1)
f4 <- fitdist(x4, "norm")
g4 <- gofstat(f4)
g4$kstest

par(mfrow=1:2)
denscomp(f3, addlegend = FALSE, main = "big sample") 
denscomp(f4, addlegend = FALSE, main = "small sample")

## ------------------------------------------------------------------------
g3$chisqtable
g3$chisqpvalue
g4$chisqtable
g4$chisqpvalue

## ------------------------------------------------------------------------
set.seed(1234)
g <- rgamma(100, shape = 2, rate = 1)
(f <- fitdist(g, "gamma"))
(f0 <- fitdist(g, "exp"))
L <- logLik(f)
k <- length(f$estimate) # number of parameters of the complete distribution
L0 <- logLik(f0)
k0 <- length(f0$estimate) # number of parameters of the simplified distribution
(stat <- 2*L - 2*L0)
(critical_value <- qchisq(0.95, df = k - k0))
(rejected <- stat > critical_value)

## ------------------------------------------------------------------------
dshiftlnorm <- function(x, mean, sigma, shift, log = FALSE) dlnorm(x+shift, mean, sigma, log=log)
pshiftlnorm <- function(q, mean, sigma, shift, log.p = FALSE) plnorm(q+shift, mean, sigma, log.p=log.p)
qshiftlnorm <- function(p, mean, sigma, shift, log.p = FALSE) qlnorm(p, mean, sigma, log.p=log.p)-shift
dshiftlnorm_no <- function(x, mean, sigma, shift) dshiftlnorm(x, mean, sigma, shift)
pshiftlnorm_no <- function(q, mean, sigma, shift) pshiftlnorm(q, mean, sigma, shift)

## ------------------------------------------------------------------------
data(dataFAQlog1)
y <- dataFAQlog1
D <- 1-min(y)
f0 <- fitdist(y+D, "lnorm")
start <- list(mean=as.numeric(f0$estimate["meanlog"]),  
              sigma=as.numeric(f0$estimate["sdlog"]), shift=D)
# works with BFGS, but not Nelder-Mead
f <- fitdist(y, "shiftlnorm", start=start, optim.method="BFGS")
summary(f)

## ---- error=FALSE--------------------------------------------------------
f2 <- try(fitdist(y, "shiftlnorm_no", start=start, optim.method="BFGS"))
print(attr(f2, "condition"))

## ------------------------------------------------------------------------
sum(log(dshiftlnorm_no(y, 0.16383978, 0.01679231, 1.17586600 )))
log(prod(dshiftlnorm_no(y, 0.16383978, 0.01679231, 1.17586600 )))
sum(dshiftlnorm(y, 0.16383978, 0.01679231, 1.17586600, TRUE ))

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  double dlnorm(double x, double meanlog, double sdlog, int give_log)
#  {
#      double y;
#  
#  #ifdef IEEE_754
#      if (ISNAN(x) || ISNAN(meanlog) || ISNAN(sdlog))
#  	return x + meanlog + sdlog;
#  #endif
#      if(sdlog <= 0) {
#  	if(sdlog < 0) ML_ERR_return_NAN;
#  	// sdlog == 0 :
#  	return (log(x) == meanlog) ? ML_POSINF : R_D__0;
#      }
#      if(x <= 0) return R_D__0;
#  
#      y = (log(x) - meanlog) / sdlog;
#      return (give_log ?
#  	    -(M_LN_SQRT_2PI   + 0.5 * y * y + log(x * sdlog)) :
#  	    M_1_SQRT_2PI * exp(-0.5 * y * y)  /	 (x * sdlog));
#      /* M_1_SQRT_2PI = 1 / sqrt(2 * pi) */
#  
#  }

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  -(M_LN_SQRT_2PI   + 0.5 * y * y + log(x * sdlog))

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  M_1_SQRT_2PI * exp(-0.5 * y * y)  /	 (x * sdlog))

## ------------------------------------------------------------------------
f2 <- fitdist(y, "shiftlnorm", start=start, lower=c(-Inf, 0, -min(y)), optim.method="Nelder-Mead")
summary(f2)
print(cbind(BFGS=f$estimate, NelderMead=f2$estimate))


## ------------------------------------------------------------------------
data(dataFAQscale1)
head(dataFAQscale1)
summary(dataFAQscale1)

## ------------------------------------------------------------------------
for(i in 6:0)
cat(10^i, try(mledist(dataFAQscale1*10^i, "cauchy")$estimate), "\n")

## ------------------------------------------------------------------------
data(dataFAQscale2)
head(dataFAQscale2)
summary(dataFAQscale2)

## ------------------------------------------------------------------------
for(i in 0:5)
cat(10^(-2*i), try(mledist(dataFAQscale2*10^(-2*i), "cauchy")$estimate), "\n")

## ----scalenormal, echo=TRUE, warning=FALSE-------------------------------
set.seed(1234)
x <- rnorm(1000, 1, 2)
fitdist(x, "norm", lower=c(-Inf, 0))

## ----shapeburr, echo=TRUE, warning=FALSE---------------------------------
x <- rburr(1000, 1, 2, 3)
fitdist(x, "burr", lower=c(0, 0, 0), start=list(shape1 = 1, shape2 = 1, 
  rate = 1))

## ----probgeom, echo=TRUE, warning=FALSE----------------------------------
x <- rgeom(1000, 1/4)
fitdist(x, "geom", lower=0, upper=1)

## ----shiftexp, echo=TRUE, warning=FALSE----------------------------------
dsexp <- function(x, rate, shift)
  dexp(x-shift, rate=rate)
psexp <- function(x, rate, shift)
  pexp(x-shift, rate=rate)
rsexp <- function(n, rate, shift)
  rexp(n, rate=rate)+shift
x <- rsexp(1000, 1/4, 1)
fitdist(x, "sexp", start=list(rate=1, shift=0), lower= c(0, -min(x)))

## ---- fig.height=3, fig.width=6------------------------------------------
pgeom(0:3, prob=1/2)
qgeom(c(0.3, 0.6, 0.9), prob=1/2)
par(mar=c(4,4,2,1), mfrow=1:2)
curve(pgeom(x, prob=1/2), 0, 10, n=301, main="c.d.f.")
curve(qgeom(x, prob=1/2), 0, 1, n=301, main="q.f.")

## ------------------------------------------------------------------------
x <- c(0, 0, 0, 0, 1, 1, 3, 2, 1, 0, 0)
median(x[-1]) #sample size 10
median(x) #sample size 11

## ---- fig.height=4, fig.width=6------------------------------------------
x <- rgeom(100, 1/3)
L2 <- function(p)
  (qgeom(1/2, p) - median(x))^2
L2(1/3) #theoretical value
curve(L2(x), 0.10, 0.95, xlab=expression(p), ylab=expression(L2(p)), main="squared differences", n=301)

## ------------------------------------------------------------------------
fitdist(x, "geom", method="qme", probs=1/2, start=list(prob=1/2), control=list(trace=1, REPORT=1))
fitdist(x, "geom", method="qme", probs=1/2, start=list(prob=1/20), control=list(trace=1, REPORT=1))

## ------------------------------------------------------------------------
fitdist(x, "geom", method="qme", probs=1/2, optim.method="SANN", start=list(prob=1/20))
fitdist(x, "geom", method="qme", probs=1/2, optim.method="SANN", start=list(prob=1/2))

## ---- fig.height=4, fig.width=6------------------------------------------
par(mar=c(4,4,2,1))
x <- rpois(100, lambda=7.5)
L2 <- function(lam)
  (qpois(1/2, lambda = lam) - median(x))^2
curve(L2(x), 6, 9, xlab=expression(lambda), ylab=expression(L2(lambda)), main="squared differences", n=201)

## ------------------------------------------------------------------------
fitdist(x, "pois", method="qme", probs=1/2, start=list(lambda=2))
fitdist(x, "pois", method="qme", probs=1/2, optim.method="SANN", start=list(lambda=2))

## ---- fig.height=4, fig.width=4, warning = FALSE-------------------------
set.seed(1234)
n <- rnorm(30, mean = 10, sd = 2)
fn <- fitdist(n, "norm")
bn <- bootdist(fn)
bn$CI
fn$estimate + cbind("estimate"= 0, "2.5%"= -1.96*fn$sd, "97.5%"= 1.96*fn$sd)
llplot(fn, back.col = FALSE)

## ---- fig.height=4, fig.width=4, warning = FALSE-------------------------
set.seed(1234)
g <- rgamma(30, shape = 0.1, rate = 10)
fg <- fitdist(g, "gamma")
bg <- bootdist(fg)
bg$CI
fg$estimate + cbind("estimate"= 0, "2.5%"= -1.96*fg$sd, "97.5%"= 1.96*fg$sd)
llplot(fg, back.col = FALSE)

## ---- fig.height=3, fig.width=4, warning = FALSE-------------------------
data(salinity)
log10LC50 <-log10(salinity)
fit <- fitdistcens(log10LC50, "norm")
# Bootstrap 
bootsample <- bootdistcens(fit, niter = 101)
#### We used only 101 iterations in that example to limit the calculation time but
#### in practice you should take at least 1001 bootstrap iterations
# Calculation of the quantile of interest (here the 5 percent hazard concentration)
(HC5 <- quantile(bootsample, probs = 0.05))
# visualizing pointwise confidence intervals on other quantiles
CIcdfplot(bootsample, CI.output = "quantile", CI.fill = "pink", xlim = c(0.5,2), main = "")

## ------------------------------------------------------------------------
exposure <- 1.2
# Bootstrap sample of the PAF at this exposure
PAF <- pnorm(exposure, mean = bootsample$estim$mean, sd = bootsample$estim$sd)
# confidence interval from 2.5 and 97.5 percentiles
quantile(PAF, probs = c(0.025, 0.975))

