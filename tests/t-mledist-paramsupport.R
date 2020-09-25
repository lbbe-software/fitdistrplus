library(fitdistrplus)



set.seed(1234)
n <- 1e2


# (1) uniform distribution fit - no fixed value

#manual check
dunif2 <- function(x, min, max)
  dunif(x, min, max)
punif2 <- function(q, min, max)
  punif(q, min, max)

x1 <- runif(n, 3, 5)

L <- function(a, b, obs)
  prod(dunif(obs, min=a, max=b))
l <- Vectorize(L, "a")
curve(l(x, b=5, obs=x1), from=1, to=3)

f1 <- fitdist(x1, "unif")
f2 <- fitdist(x1, "unif2", start=list(min=0, max=10), lower=c(-Inf, max(x1)),
              upper=c(min(x1), Inf))
c(logLik(f1), logLik(f2))

delta <- .2
llsurface(x1, "unif", plot.arg = c("min", "max"), min.arg=c(1, 5-delta),
          max.arg=c(3+delta, 7), main="likelihood surface for uniform")
abline(v=3, h=5, col="red", lty=2)
points(f1$estimate[1], f1$estimate[2], pch="x", col="violet")
points(f2$estimate[1], f2$estimate[2], pch="o", col="blue")

# (2) uniform distribution fit - fixed value

f3 <- fitdist(x1, "unif", fix.arg=list(min=2.5))
logLik(f3)

llcurve(x1, "unif", plot.arg="max", min.arg = 5-delta, max.arg=7)

f4 <- fitdist(x1, "unif", fix.arg=list(max=5.5))
logLik(f4)


# (3) four parameter beta - also known as PERT distribution
require(mc2d)

x2 <- rpert(n, 0, 1, 2, 3)

f1 <- fitdist(x2, "pert", start=list(min=-1, mode=0, max=10, shape=1),
              lower=c(-Inf, -Inf, max(x2), 0),
              upper=c(min(x2), Inf, Inf, Inf))
f2 <- fitdist(x2, "pert", start=list(mode=1, shape=1), 
              lower=c(-Inf, 0),
              upper=c(Inf, Inf),
              fix.arg=list(min=min(x2)-1e-6, max=max(x2)+1e-6))

gofstat(list(f1,f2))
cdfcomp(list(f1,f2))
