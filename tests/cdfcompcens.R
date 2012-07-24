library(fitdistrplus)

# (1) Plot various distributions fitted to bacterial contamination data
#
data(smokedfish)
fitsfn <- fitdistcens(smokedfish,"norm")
summary(fitsfn)

fitsfl <- fitdistcens(smokedfish,"logis")
summary(fitsfl)

dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q,a,b) exp(-exp((a-q)/b))
qgumbel <- function(p,a,b) a-b*log(-log(p))
fitsfg<-fitdistcens(smokedfish,"gumbel",start=list(a=-3,b=3))
summary(fitsfg)

cdfcompcens(list(fitsfn,fitsfl,fitsfg))
cdfcompcens(list(fitsfn,fitsfl,fitsfg),datacol="orange",
legendtext=c("normal","logistic","Gumbel"),
main="bacterial contamination fits",
xlab="bacterial concentration (CFU/g)",ylab="F")
cdfcompcens(list(fitsfn,fitsfl,fitsfg),datacol="orange",
legendtext=c("normal","logistic","Gumbel"),
main="bacterial contamination fits",
xlab="bacterial concentration (CFU/g)",ylab="F",
xlegend = "center")

# (2) Same plot in x logscale
#
cdfcompcens(list(fitsfn,fitsfl,fitsfg),xlog = TRUE)
