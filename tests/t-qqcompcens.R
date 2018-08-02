library(fitdistrplus)

data(smokedfish)
fitsf  <-  fitdistcens(smokedfish,"lnorm")
plot(fitsf)
qqcompcens(fitsf)
qqcompcens(fitsf, fillrect = NA)
qqcompcens(fitsf, fitcol = "black")
qqcompcens(fitsf, fitcol = "black", fillrect = NA)
qqcompcens(fitsf, ylim = c(0,150)) 
qqcompcens(fitsf, xlim = c(0,150)) 
qqcompcens(fitsf, xlim = c(0,150), ylim = c(0, 120)) 

if (requireNamespace("ggplot2", quietly = TRUE)) {
  qqcompcens(fitsf, plotstyle = "ggplot")
  qqcompcens(fitsf, fillrect = NA, plotstyle = "ggplot")
  qqcompcens(fitsf, fitcol = "black", plotstyle = "ggplot")
  qqcompcens(fitsf, fitcol = "black", fillrect = NA, plotstyle = "ggplot")
  qqcompcens(fitsf, ylim = c(0,150), plotstyle = "ggplot")
  qqcompcens(fitsf, xlim = c(0,150), plotstyle = "ggplot")
  qqcompcens(fitsf, xlim = c(0,150), ylim = c(0, 120), plotstyle = "ggplot")
}

data(fluazinam)
log10EC50 <-log10(fluazinam)
fln <- fitdistcens(log10EC50,"norm")
plot(fln)
qqcompcens(fln)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  qqcompcens(fln, plotstyle = "ggplot")
}

data(salinity)
log10LC50 <-log10(salinity)
plotdistcens(log10LC50)
plotdistcens(log10LC50, NPMLE = FALSE)
fn <- fitdistcens(log10LC50,"norm")
fl <- fitdistcens(log10LC50,"logis")
plot(fn)
plot(fl)
qqcompcens(fn)
qqcompcens(fl)
qqcompcens(list(fn, fl))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  qqcompcens(fn, plotstyle = "ggplot")
  qqcompcens(fl, plotstyle = "ggplot")
  qqcompcens(list(fn, fl), plotstyle = "ggplot")
}


require(actuar)
data(salinity)
fln <- fitdistcens(salinity,"lnorm")
fll <- fitdistcens(salinity,"llogis")
plot(fln)
par(mfrow = c(2,1))
qqcompcens(fln)
qqcompcens(fll)
par(mfrow = c(1,1))
qqcompcens(list(fln, fll))
qqcompcens(list(fln, fll), ynoise = FALSE)
qqcompcens(list(fln, fll), fitcol = c("blue", "orange"))
qqcompcens(list(fln, fll), xlogscale = TRUE, ylogscale = TRUE)
qqcompcens(list(fln, fll), ylogscale = TRUE) 
qqcompcens(list(fln, fll), xlogscale = TRUE, ynoise = FALSE) 

if (requireNamespace("ggplot2", quietly = TRUE)) {
  qqcompcens(list(fln, fll), plotstyle = "ggplot")
  qqcompcens(list(fln, fll), ynoise = FALSE, plotstyle = "ggplot")
  qqcompcens(list(fln, fll), fitcol = c("blue", "orange"), plotstyle = "ggplot")
  qqcompcens(list(fln, fll), xlogscale = TRUE, ylogscale = TRUE, plotstyle = "ggplot")
  qqcompcens(list(fln, fll), ylogscale = TRUE, plotstyle = "ggplot")
  qqcompcens(list(fln, fll), xlogscale = TRUE, ynoise = FALSE, plotstyle = "ggplot")
}
