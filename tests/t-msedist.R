library(fitdistrplus)
nboot <- 1000
nboot <- 10

set.seed(123)

#manual implementation
geomeanspacing <- function(pdist, obs, ...)
{
  sx <- c(-Inf, sort(obs), Inf)
  n <- length(sx)
  Di <- pdist(sx[-1], ...) - pdist(sx[-n], ...)
  mean(log(Di))
}

#--------------------------------------------------------
# exponential sample

x1 <- rexp(1e3)
lam <- seq(0.1, 2, length=51)
Sn <- sapply(lam, function(x) geomeanspacing(pexp, obs=x1, rate=x))

Dn <- function(theta)
  -geomeanspacing(pexp, obs=x1, rate=theta[1])
#theoretical optimum
Dn(1/mean(x1))

#check curve behavior
par(mar=c(4,4,2,1))
plot(lam, Sn, type="l", xlab="theta", main="Average spacing logarithm")
abline(v=1, col="green")
abline(v=msedist(x1, "exp")$estimate, lty=2, col="blue")
legend("bottomright", lty=1:2, col=c("green", "blue"), leg=c("theoretical value", "fitted value"))

msedist(x1, "exp", control=list(trace=0, REPORT=1))

mse_exp <- fitdist(x1, "exp", method="mse")
plot(mse_exp)
summary(mse_exp)
gofstat(mse_exp)

mse_exp_boot <- bootdist(mse_exp, niter = nboot)
plot(mse_exp_boot)
abline(v=1, col="green")
abline(v=msedist(x1, "exp")$estimate, lty=2, col="blue")
legend("bottomright", lty=1:2, col=c("green", "blue"), leg=c("theoretical value", "fitted value"))


# library(BMT)
# x <- rBMT(1e3, 1/4, 3/4)
# 
# BMTfit.mpse(x)
# fitdist(x, "BMT", method="mse", start=list(p3=1/2, p4=1/2, p1=-1/2, p2=1/2), lower=c(0, 0, -Inf, 0), 
#         upper=c(1,1,0, Inf))
# pBMT(x, p3=1/2, p4=1/2, p1=-1/2, p2=1/2)

#--------------------------------------------------------
# lognormal sample


x1 <- rlnorm(1e3, 0, 1)
mu <- seq(-1, 1, length=51)
Sn <- sapply(mu, 
             function(x) geomeanspacing(plnorm, obs=x1, mean=x, sd=1))


Dn <- function(theta)
  -geomeanspacing(plnorm, obs=x1, mean=theta[1], sd=theta[2])

plot(mu, Sn, type="l")
abline(v=0)


optim(c(2,2), Dn)
msedist(x1, "lnorm", control=list(trace=0, REPORT=1))


mse_lnorm <- fitdist(x1, "lnorm", method="mse")
mle_lnorm <- fitdist(x1, "lnorm", method="mle")
plot(mse_lnorm)
summary(mse_lnorm)
cdfcomp(list(mse_lnorm, mle_lnorm))
gofstat(list(mse_lnorm, mle_lnorm))
mse_lnorm_boot <- bootdist(mse_lnorm, niter = nboot)
par(mar=c(4,4,2,1))
plot(mse_lnorm_boot, enhance = TRUE, trueval=c(0,1))

#--------------------------------------------------------
# Pareto sample

library(actuar)

x1 <- rburr(1e4, 2,2,2)

Dn <- function(theta)
  -geomeanspacing(pburr, obs=x1, shape1=theta[1], shape2=theta[2], rate=theta[3])
Dn(c(1,1,10))

optim(c(1,1,10), Dn)

msedist(x1, "burr", start=list(shape1=1, shape2=1, rate=10), control=list(trace=0, REPORT=1))


mse_burr <- fitdist(x1, "burr", method="mse", start=list(shape1=1, shape2=1, rate=10))
mle_burr <- fitdist(x1, "burr", method="mle", start=list(shape1=1, shape2=1, rate=10))
plot(mse_burr)
summary(mse_burr)
cdfcomp(list(mse_burr, mle_burr))
gofstat(list(mse_burr, mle_burr))
mse_burr_boot <- bootdist(mse_burr, niter = pmin(nboot,100))
plot(mse_burr_boot, enhance = TRUE, trueval=c(2,2,2))



#--------------------------------------------------------
# Poisson sample

x1 <- rpois(1e3, 15)

geomeanSpacingUnique <- function(pdist, obs, ...)
{
  sx <- c(-Inf, unique(sort(obs)), Inf)
  n <- length(sx)
  Di <- pdist(sx[-1], ...) - pdist(sx[-n], ...)
  mean(log(Di))
}

geomeanSpacingWeight <- function(pdist, obs, weights, ...)
{
  sx <- c(-Inf, unique(sort(obs)), Inf)
  weights <- c(1, weights)
  n <- length(sx)
  Di <- pdist(sx[-1], ...) - pdist(sx[-n], ...)
  mean(weights*log(Di))
}

DnUnique <- function(theta)
  -geomeanSpacingUnique(ppois, obs=x1, lambda=theta[1])

DnWeight <- function(theta, weights)
  -geomeanSpacingWeight(ppois, obs=x1, lambda=theta[1], weights=weights)


optimize(DnWeight, c(1, 30), weights=as.numeric(table(x1)))
optimize(DnUnique, c(1, 30))
optimize(Dn, c(1, 30)) #does not converge

mle_pois1 <- fitdist(x1, "pois", method="mle")
#no weight
mse_pois1 <- fitdist(x1, "pois", method="mse")
#with weight
mse_pois2 <- fitdist(unique(sort(x1)), "pois", method="mse", weights=as.numeric(table(x1)))
plot(mse_pois1)
plot(mse_pois2)
summary(mse_pois1)
gofstat(mse_pois1)
gofstat(mse_pois2)

par(mfrow=c(1,1))
cdfcomp(list(mle_pois1, mse_pois1), addlegend = FALSE, fitlty = 1)
curve(ppois(x, lambda=mse_pois2$estimate), type="s", col="blue", add=TRUE)
legend("bottomright", lty=1, col=c("red", "green", "blue"), leg=c("MLE", "MSE no weight", "MSE with weight"))



#--------------------------------------------------------
# real dataset
# library(CASdatasets)
# data("ushustormloss")
# x <- ushustormloss$Normalized.CL05
# 
# plot(Normalized.CL05 ~ Year, data=ushustormloss, type="h", main="Normalized Hurricane Damages in United States")
# 
# mse_burr <- fitdist(x, "burr", method="mse", start=list(shape1=1, shape2=1, rate=10), lower=0)
# mle_burr0 <- fitdist(x, "burr", method="mle", start=list(shape1=1, shape2=1, rate=10), lower=0)
# 
# cbind(MSE=coef(mse_burr), MLE=coef(mle_burr0))
# 
# 
# setwd("~/Desktop/")
# par(mar=c(4,4,2,1))
# pdf("Ushustorm-cdfcomp.pdf", 6, 6)
# cdfcomp(list(mse_burr, mle_burr0), xlogscale = TRUE, do.points = FALSE)
# dev.off()
# pdf("Ushustorm-qqcomp.pdf", 6, 6)
# qqcomp(list(mse_burr, mle_burr0), xlogscale=TRUE, ylogscale=TRUE)
# dev.off()
# pdf("Ushustorm-ppcomp.pdf", 6, 6)
# ppcomp(list(mse_burr, mle_burr0))
# dev.off()
# 
# gofstat(list(mse_burr, mle_burr0))
# 
# mse_iburr <- fitdist(x, "invburr", method="mse", start=list(shape1=1, shape2=1, rate=10), lower=0)
# mle_iburr0 <- fitdist(x, "invburr", method="mle", start=list(shape1=1, shape2=1, rate=10), lower=0)
# 
# gofstat(list(mse_iburr, mle_iburr0))
# cdfcomp(list(mse_iburr, mle_iburr0))

