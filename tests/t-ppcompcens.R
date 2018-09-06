library(fitdistrplus)

visualize <- FALSE # TRUE for manual tests with visualization of results

data(smokedfish)
fitsf  <-  fitdistcens(smokedfish,"lnorm")
plot(fitsf)
ppcompcens(fitsf)
ppcompcens(fitsf, fillrect = NA)
ppcompcens(fitsf, fitcol = "black")
ppcompcens(fitsf, fitcol = "black", fillrect = NA)
ppcompcens(fitsf, ylim = c(0.4,1)) 
ppcompcens(fitsf, xlim = c(0.4,1)) 
ppcompcens(fitsf, xlim = c(0.4,1), ylim = c(0,1)) 
ppcompcens(fitsf, xlim = c(0.5,0.99), xlogscale = TRUE)
try(ppcompcens(fitsf, xlogscale = TRUE))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ppcompcens(fitsf, plotstyle = "ggplot")
}
if (requireNamespace("ggplot2", quietly = TRUE) & visualize) {
  ppcompcens(fitsf, fillrect = NA, plotstyle = "ggplot")
  ppcompcens(fitsf, fitcol = "black", plotstyle = "ggplot")
  ppcompcens(fitsf, fitcol = "black", fillrect = NA, plotstyle = "ggplot")
  ppcompcens(fitsf, ylim = c(0.4,1), plotstyle = "ggplot")
  ppcompcens(fitsf, xlim = c(0.4,1), plotstyle = "ggplot")
  ppcompcens(fitsf, xlim = c(0.4,1), ylim = c(0,1), plotstyle = "ggplot")
  ppcompcens(fitsf, xlim = c(0.5,0.99), xlogscale = TRUE, plotstyle = "ggplot")
}

if (visualize)
{
  data(fluazinam)
  log10EC50 <-log10(fluazinam)
  fln <- fitdistcens(log10EC50,"norm")
  plot(fln)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ppcompcens(fln, plotstyle = "ggplot")
  }
}

data(salinity)
log10LC50 <-log10(salinity)
fn <- fitdistcens(log10LC50,"norm")
fl <- fitdistcens(log10LC50,"logis")
ppcompcens(list(fn, fl))

if (visualize)
{
  plotdistcens(log10LC50)
  plotdistcens(log10LC50, NPMLE = FALSE)
  plot(fn)
  plot(fl)
  ppcompcens(fn)
  ppcompcens(fl)
  ppcompcens(list(fn, fl), ynoise = FALSE)
  ppcompcens(list(fn, fl), xlogscale = TRUE, xlim = c(0.01, 0.6))
}

if (requireNamespace ("ggplot2", quietly = TRUE) ) {
  ppcompcens(list(fn, fl), plotstyle = "ggplot", fitcol = "red")
}


if (requireNamespace ("ggplot2", quietly = TRUE) & visualize) {
  ppcompcens(fn, plotstyle = "ggplot")
  ppcompcens(fl, plotstyle = "ggplot")
  ppcompcens(list(fn, fl), plotstyle = "ggplot", fitcol = "red")
  ppcompcens(list(fn, fl), ynoise = FALSE, plotstyle = "ggplot")
  ppcompcens(list(fn, fl), xlogscale = TRUE, xlim = c(0.01, 0.6), plotstyle = "ggplot")
}

