library(actuar)

n <- 1e3
x <- rztpois(n, 2)

hist(x)
hist(rpois(n, 2)+1)

cdfcomp(fitdist(x-1, "pois"))
