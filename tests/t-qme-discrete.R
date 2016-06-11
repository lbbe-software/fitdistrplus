library(fitdistrplus)
set.seed(1234)

#Poisson

x <- rpois(100, lambda=7.5)
L2 <- function(lam)
  (qpois(1/2, lambda = lam) - median(x))^2
curve(L2(x), 5, 9, xlab=expression(lambda), ylab=expression(L2(lambda)), main="squared differences", n=201)

fitdist(x, "pois", method="qme", probs=1/2, start=list(lambda=2), control=list(trace=1, REPORT=1))
fitdist(x, "pois", method="qme", probs=1/2, start=list(lambda=6.8), control=list(trace=1, REPORT=1))
fitdist(x, "pois", method="qme", probs=1/2, start=list(lambda=15), control=list(trace=1, REPORT=1))

fitdist(x, "pois", method="qme", optim.method="SANN", probs=1/2, start=list(lambda=2), control=list(trace=1, REPORT=100))
fitdist(x, "pois", method="qme", optim.method="SANN", probs=1/2, start=list(lambda=17), control=list(trace=1, REPORT=100))

if(FALSE) #does not work
{
# library(nleqslv)
# 
# F <- function(x, obs)
#   qpois(1/2, lambda = x) - median(obs)
# nleqslv(2, F, obs=x, control=list(allowSingular=TRUE))
# nleqslv(6.8, F, obs=x, control=list(allowSingular=TRUE))

}


#Geometric

x <- rgeom(100, 1/3)
L2 <- function(p)
  (qgeom(1/2, p) - median(x))^2
curve(L2(x), 0.10, 0.95, xlab=expression(p), ylab=expression(L2(p)), main="squared differences", n=301)

fitdist(x, "geom", method="qme", probs=1/2, start=list(prob=1/2), control=list(trace=1, REPORT=1))
fitdist(x, "geom", method="qme", probs=1/2, start=list(prob=1/20), control=list(trace=1, REPORT=1))
