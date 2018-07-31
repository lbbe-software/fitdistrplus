library(fitdistrplus)

# (1) Plot various distributions fitted to bacterial contamination data
#
data(smokedfish)
Clog10 <- log10(smokedfish)

fitsfn <- fitdistcens(Clog10,"norm")
fitsfl <- fitdistcens(Clog10,"logis")

dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q,a,b) exp(-exp((a-q)/b))
qgumbel <- function(p,a,b) a-b*log(-log(p))
fitsfg <- fitdistcens(Clog10, "gumbel", start=list(a=-3,b=3))

cdfcompcens(list(fitsfn,fitsfl,fitsfg))
cdfcompcens(list(fitsfn,fitsfl,fitsfg), datacol="grey",
    legendtext=c("normal","logistic","Gumbel"),
    main="bacterial contamination fits",
    xlab="bacterial concentration (CFU/g)", ylab="F",
    xlegend = "center", lines01 = TRUE)

cdfcompcens(list(fitsfn, fitsfl, fitsfg), NPMLE.method = "Turnbull")

# Same plot in y logscale
cdfcompcens(list(fitsfn, fitsfl, fitsfg), NPMLE.method = "Turnbull",
            ylogscale = TRUE, ylim=c(.5, .99))
cdfcompcens(list(fitsfn, fitsfl, fitsfg), NPMLE.method = "Wang",
            ylogscale = TRUE, ylim=c(.5, .99))

if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcompcens(list(fitsfn,fitsfl,fitsfg), plotstyle = "ggplot")
  cdfcompcens(list(fitsfn,fitsfl,fitsfg), datacol="grey",
              legendtext=c("normal","logistic","Gumbel"),
              main="bacterial contamination fits",
              xlab="bacterial concentration (CFU/g)", ylab="F",
              xlegend = "center", lines01 = TRUE, plotstyle = "ggplot")
  
  # Same plot in y logscale
  cdfcompcens(list(fitsfn, fitsfl, fitsfg), NPMLE.method = "Wang",
              ylogscale = TRUE, ylim=c(.5, .99), plotstyle = "ggplot")
}

# Use of x logscale
if(any(installed.packages()[,"Package"] == "actuar"))
{
  require(actuar)
  data(smokedfish)
  fln <- fitdistcens(smokedfish,"lnorm")
  fll <- fitdistcens(smokedfish,"llogis")
  cdfcompcens(list(fln, fll))
  cdfcompcens(list(fln, fll), xlogscale = TRUE)
  cdfcompcens(list(fln, fll), xlogscale = TRUE, xlim = c(0.01, 1000))
  cdfcompcens(list(fln, fll), NPMLE.method = "Turnbull",
              xlogscale = TRUE, xlim = c(0.01, 1000))
  
  if (requireNamespace ("ggplot2", quietly = TRUE)) {
    cdfcompcens(list(fln, fll), plotstyle = "ggplot")
    cdfcompcens(list(fln, fll), xlogscale = TRUE, plotstyle = "ggplot")
    cdfcompcens(list(fln, fll), xlogscale = TRUE, xlim = c(0.01, 1000), plotstyle = "ggplot")
  }
}
# same plot using argument add
cdfcompcens(fitsfn, addlegend = FALSE, fitcol = "red")
cdfcompcens(fitsfl, addlegend = FALSE, fillrect = NA, fitcol = "green", add = TRUE)
cdfcompcens(fitsfg, addlegend = FALSE, fillrect = NA, fitcol = "blue", add = TRUE)

cdfcompcens(list(fitsfn, fitsfl, fitsfg), addlegend = FALSE, fitcol = 2:4, fitlty = 1, plotstyle = "ggplot")
