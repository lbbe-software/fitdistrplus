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
