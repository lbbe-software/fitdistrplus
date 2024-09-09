require(fitdistrplus)
data(salinity)
log10LC50 <-log10(salinity)
fitn <- fitdistcens(log10LC50, "norm")
fitc <- fitdistcens(log10LC50, "cauchy")
cdfcompcens(list(fitc, fitn), xlim = c(1, 1.8), main = "",
            fitcol = c("blue", "orange"), 
            fitlty = c(1, 1), fitlwd = c(2, 2))

data(salinity)
fitln <- fitdistcens(salinity, "lnorm")
fitw <- fitdistcens(salinity, "weibull")
fitg <- fitdistcens(salinity, "gamma")
cdfcompcens(list(fitln, fitw, fitg), main = "",
            fitcol = c("blue", "orange", "green"), 
            fitlty = c(1, 1, 1), fitlwd = c(2, 2, 2))
cdfcompcens(list(fitln, fitw, fitg), main = "",
            fitcol = c("blue", "orange", "green"), 
            fitlty = c(1, 1, 1), fitlwd = c(2, 2, 2),
            xlogscale = TRUE)
