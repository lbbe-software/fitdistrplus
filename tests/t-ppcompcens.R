library(fitdistrplus)

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

data(fluazinam)
log10EC50 <-log10(fluazinam)
fln <- fitdistcens(log10EC50,"norm")
plot(fln)
ppcompcens(fln)


data(salinity)
log10LC50 <-log10(salinity)
plotdistcens(log10LC50)
plotdistcens(log10LC50, NPMLE = FALSE)
fn <- fitdistcens(log10LC50,"norm")
fl <- fitdistcens(log10LC50,"logis")
plot(fn)
plot(fl)
ppcompcens(fn)
ppcompcens(fl)
ppcompcens(list(fn, fl))
ppcompcens(list(fn, fl), ynoise = FALSE)
ppcompcens(list(fn, fl), xlogscale = TRUE, xlim = c(0.01, 0.6))



