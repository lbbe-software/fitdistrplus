library(fitdistrplus)

nbboot <- 1001
nbboot <- 10

# (1) Fit of a gamma distribution to serving size data
# using default method (maximum likelihood estimation)
# followed by parametric bootstrap
#

s1 <- rgamma(50, 3, 2)
f1 <- fitdist(s1, "gamma")
b1 <- bootdist(f1, niter=nbboot, silent=TRUE)

plot(b1)
quantile(b1)

str(b1)
b1$estim

b <- b1
qdistname <- paste("q", b$fitpart$distname, sep="")
pdistname <- paste("p", b$fitpart$distname, sep="")

Fp <- function(x)
{  
calcp <- function(i)
{
  parai <- c(as.list(b$estim[i, ]), as.list(b$fitpart$fix.arg))
  do.call(pdistname, c(list(q=x), as.list(parai)))
}
res <- t(sapply(1:b$nbboot, calcp))
rownames(res) <- 1:b$nbboot
colnames(res) <- paste0("x=", x)
res
}

x <- seq(0, max(b$fitpart$data), length=101)

cdfcomp(b1$fitpart)


matlines(x, CIband, col="red", lty=2)

par(mfrow=c(1,3))
cdfband(b1, CI.level=90/100)
cdfband(b1, CI.level=90/100, CI.type = "less")
cdfband(b1, CI.level=90/100, CI.type = "greater")


par(mfrow=c(1,3))
cdfband(b1, CI.level=90/100, CI.col="grey90", CI.fill=TRUE, datacol="blue", verticals=TRUE)
cdfband(b1, CI.level=90/100, CI.col="grey90", CI.type = "less", CI.fill=TRUE, datacol="blue", do.points=FALSE)
cdfband(b1, CI.level=90/100, CI.col="grey90", CI.type = "greater", CI.fill=TRUE, datacol="blue", verticals=TRUE)


x <- seq(min(b$fitpart$data)-1, max(b$fitpart$data)+1, length=101)

CIband <- apply(Fp(x), 2, quantile, probs=.99)
plot(ecdf(s1))
polygon(c(x, max(x)*1.2, max(x)*1.2), c(CIband, 1, 0), col="grey90")
plot(ecdf(s1), add=TRUE)

CIband <- apply(Fp(x), 2, quantile, probs=1-.99)
plot(ecdf(s1))
polygon(c(x, min(x)-1, min(x)-1), c(CIband, 1, 0), col="grey90")
plot(ecdf(s1), add=TRUE)


CIband <- t(apply(Fp(x), 2, quantile, probs=c(.025, .975)))
plot(ecdf(s1))
polygon(c(x, rev(x)), c(CIband[,2], rev(CIband[,1])), col="grey90")
plot(ecdf(s1), add=TRUE)


