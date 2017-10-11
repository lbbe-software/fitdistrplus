library(fitdistrplus)

# (1) Plot various distributions fitted to bacterial contamination data
#
data(smokedfish)
Clog10 <- log10(smokedfish)

fitsfn <- fitdistcens(Clog10,"norm")
summary(fitsfn)

fitsfl <- fitdistcens(Clog10,"logis")
summary(fitsfl)

dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q,a,b) exp(-exp((a-q)/b))
qgumbel <- function(p,a,b) a-b*log(-log(p))
fitsfg<-fitdistcens(Clog10, "gumbel", start=list(a=-3,b=3))
summary(fitsfg)

cdfcompcens(list(fitsfn,fitsfl,fitsfg))
cdfcompcens(list(fitsfn,fitsfl,fitsfg), datacol="orange",
    legendtext=c("normal","logistic","Gumbel"),
    main="bacterial contamination fits",
    xlab="bacterial concentration (CFU/g)", ylab="F",
    xlegend = "center", lines01 = TRUE)

# Same plot in y logscale
cdfcompcens(list(fitsfn, fitsfl, fitsfg), ylog = TRUE)

# same plot using argument add
cdfcompcens(fitsfn, addlegend = FALSE, fitcol = "red")
cdfcompcens(fitsfl, addlegend = FALSE, fitcol = "green", add = TRUE)
cdfcompcens(fitsfg, addlegend = FALSE, fitcol = "blue", add = TRUE)


