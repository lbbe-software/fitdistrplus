library(fitdistrplus)




# (1) basic fit of a normal distribution with maximum likelihood estimation
#

set.seed(1234)
x1 <- rnorm(n=100)

#correct usage
fitdistrplus:::mledist2(x1,"norm")
fitdistrplus:::mledist2(x1,"norm", start=function(x) list(mean=0, sd=1))
fitdistrplus:::mledist2(x1,"norm", fix.arg=list(mean=1/2))
fitdistrplus:::mledist2(x1,"norm", fix.arg=list(mean=1/2), start=list(sd=1))
fitdistrplus:::mledist2(x1,"norm", fix.arg=function(x) list(mean=0), start=list(sd=1))
fitdistrplus:::mledist2(x1,"norm", fix.arg=function(x) list(mean=mean(x)), start=list(sd=1))

#wrong usage
try( fitdistrplus:::mledist2(x1,"norm", start=list(a=1/2)) )
try( fitdistrplus:::mledist2(x1,"norm", start=function(x) list(a=0, b=1)) )
try( fitdistrplus:::mledist2(x1,"norm", fix.arg=list(a=1/2)) )
try( fitdistrplus:::mledist2(x1,"norm", fix.arg=function(x) list(a=0), start=list(sd=1)) )


#relevant example
#zero modified geometric distribution
dzmgeom <- function(x, p1, p2)
{
  p1 * (x == 0) + (1-p1)*dgeom(x-1, p2)
}
rzmgeom <- function(n, p1, p2)
{
  u <- rbinom(n, 1, 1-p1) #prob to get zero is p1
  print(table(u))
  u[u != 0] <- rgeom(sum(u != 0), p2)+1
  u
}
# check
# dzmgeom(0:5, 1/2, 1/10)



x2 <- rzmgeom(1000, 1/2, 1/10)
x3 <- rzmgeom(1000, 1/3, 1/10)
x4 <- rzmgeom(1000, 1/4, 1/10)
table(x2)
#this is the MLE which converges almost surely and in distribution.
initp1 <- function(x) list(p1=mean(x == 0))

fitdistrplus:::mledist2(x2, "zmgeom", fix.arg=initp1, start=list(p2=1/2))[c("estimate", "fix.arg")]
fitdistrplus:::mledist2(x3, "zmgeom", fix.arg=initp1, start=list(p2=1/2))[c("estimate", "fix.arg")]
fitdistrplus:::mledist2(x4, "zmgeom", fix.arg=initp1, start=list(p2=1/2))[c("estimate", "fix.arg")]

