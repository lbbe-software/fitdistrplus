require(fitdistrplus)

## first example
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


## second example
data(groundbeef)
x1 <- groundbeef$serving[1:50]
f1 <- fitdist(x1, "gamma")
b1 <- bootdist(f1, niter=51)
# print(b1)
# plot(b1)
# plot(b1, enhance=TRUE)
# summary(b1)
# quantile(b1)
par(bg = "black")
CIcdfplot(b1, CI.output = "quantile", xlim = c(10, 140),
          CI.col = "red",fitlwd = 2, datacol = "yellow", CI.fill = "orange")

CIcdfplot(b1, CI.output = "quantile", xlim = c(10, 140),
          CI.col = "blue",fitlwd = 2, datacol = "snow",
          fitcol = "blue", CI.fill = "lightblue")


data(groundbeef)
x1 <- groundbeef$serving[1:30]
f1 <- fitdist(x1, "gamma")
b1 <- bootdist(f1, niter=51)
# print(b1)
# plot(b1)
# plot(b1, enhance=TRUE)
# summary(b1)
# quantile(b1)
par(bg = "black")
CIcdfplot(b1, CI.output = "quantile", xlim = c(-10, 140),
          CI.col = "red",fitlwd = 2, datacol = "yellow", CI.fill = "orange")

