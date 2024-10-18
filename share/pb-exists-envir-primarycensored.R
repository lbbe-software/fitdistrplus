
install.packages(
  "primarycensored",
  repos = "https://epinowcast.r-universe.dev"
)

library(primarycensored)
?fitdistdoublecens

set.seed(123)
n <- 1000
true_mean <- 5
true_sd <- 2
pwindow <- 2
swindow <- 2
D <- 10
samples <- rprimarycensored(
  n, rnorm,
  mean = true_mean, sd = true_sd,
  pwindow = pwindow, swindow = swindow, D = D
)
#etrange comme tirage
table(samples)

delay_data <- data.frame(
  left = samples,
  right = samples + swindow
)

fit_norm <- fitdistdoublecens(
  delay_data,
  distr = "norm",
  start = list(mean = 0, sd = 1),
  D = D, pwindow = pwindow
)

head(delay_data)

curve(pprimarycensored(x, pnorm, mean = 0, sd = 1, D=D, pwindow=pwindow), from=-5, to=10)
lines(ecdf(samples))

#c'est tres voisin de fitdistcens...
fitdistcens(delay_data, "norm")

fitdistcens(delay_data, dnorm)


