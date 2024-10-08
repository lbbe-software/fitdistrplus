require("fitdistrplus")
nsample <- 10

# (1) Plot of an empirical distribution with changing 
# of default line types for CDF and colors
# and optionally adding a density line
#
set.seed(1234)
x1 <- rnorm(n=nsample)
plotdist(x1)
plotdist(x1, col="blue", type="b", pch=16)
plotdist(x1, type="s")
plotdist(x1, demp = TRUE)
plotdist(x1, demp = TRUE, adjust = 1.5)
plotdist(x1, breaks = 2)
plotdist(x1, breaks = seq(-3, 2))

# (2) Plot of a discrete distribution against data
#
set.seed(1234)
x2 <- rpois(n=nsample, lambda = 2)
plotdist(x2, discrete=TRUE)
plotdist(x2, "pois", para=list(lambda = mean(x2)))
plotdist(x2, "pois", para=list(lambda = mean(x2)), lwd="2")

# (3) Plot of a continuous distribution against data
#
xn <- rnorm(n=nsample, mean=10, sd=5)
plotdist(xn, "norm", para=list(mean=mean(xn), sd=sd(xn)))
plotdist(xn, "norm", para=list(mean=mean(xn), sd=sd(xn)), pch=16)
plotdist(xn, "norm", para=list(mean=mean(xn), sd=sd(xn)), pch=16,
         demp = TRUE)

# (4) Plot of serving size data
#
data(groundbeef)
plotdist(groundbeef$serving, type="s")

# (5) Plot of numbers of parasites with a Poisson distribution 
data(toxocara)
number <- toxocara$number
plotdist(number, discrete = TRUE)
plotdist(number,"pois",para=list(lambda=mean(number)))

