require(fitdistrplus)
data(salinity)
log10LC50 <-log10(salinity)
fit <- fitdistcens(log10LC50, "norm")
