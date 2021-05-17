library(fitdistrplus)

compare.plotdistcens <- function(d)
{
  par(mfrow = c(2,2))
  plotdistcens(d, NPMLE.method = "Turnbull.middlepoints")
  plotdistcens(d, NPMLE.method = "Turnbull.intervals")
  plotdistcens(d, NPMLE.method = "Wang")
  
}

compare.npmle <- function(d)
{
  npmleT <- fitdistrplus:::npmle(d, method = "Turnbull.intervals")
  npmleW <- fitdistrplus:::npmle(d, method = "Wang")
  print(npmleT)
  print(npmleW)
  cat("nb of intervals for Turnbull: ", nrow(npmleT), "\n")
  cat("nb of intervals for Wang: ", nrow(npmleW), "\n")
  xmin <- min(c(npmleT$left[is.finite(npmleT$left)],npmleT$left[is.finite(npmleT$left)]))
  xmax <- max(c(npmleT$right[is.finite(npmleT$right)],npmleT$right[is.finite(npmleT$right)]))
  par(mfrow = c(2, 1))
  plot(npmleT$left, npmleT$p, xlim = c(xmin, xmax), main = "left")
  points(npmleW$left, npmleW$p, col = "red", pch = 4)
  plot(npmleT$right, npmleT$p, xlim = c(xmin, xmax), main = "right")
  points(npmleW$right, npmleW$p, col = "red", pch = 4)
}
#### Comparison of plotdistcens with different NPMLE methods
if(FALSE)
{
  
  # d1 = trivial case with only interval censored data
  d1 <- data.frame(left = c(1, 2, 3, 4, 3, 7), right = c(2, 5, 3, 7, 8, 9))
  d <- d1
  par(mfrow = c(2,2))
  compare.plotdistcens(d)
  compare.npmle(d)
  fa <- fitdistcens(d, "norm")
  fb <- fitdistcens(d, "logis")
  par(mfrow = c(2,2))
  cdfcompcens(list(fa, fb), NPMLE.method = "Turnbull.middlepoints")
  cdfcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  cdfcompcens(list(fa, fb), NPMLE.method = "Wang")
  par(mfrow = c(1,2))
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  
  
  # d2 = case with left and right censored data
  data(smokedfish)
  d2 <- smokedfish
  d <- d2
  compare.plotdistcens(d)
  compare.npmle(d)
  fa <- fitdistcens(d, "lnorm")
  fb <- fitdistcens(d, "gamma")
  par(mfrow = c(2,2))
  cdfcompcens(list(fa, fb), NPMLE.method = "Turnbull.middlepoints")
  cdfcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  cdfcompcens(list(fa, fb), NPMLE.method = "Wang")
  par(mfrow = c(1,2))
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  
  # d3 = case with also rigth censored data
  d3 <- data.frame(left = c(-1.4, 1.18, -1.4, 2, -1.4, 0),
                   right = c(1, 1.18, 2, NA, 0, 2))
  d <- d3
  compare.plotdistcens(d)
  compare.npmle(d)
  fa <- fitdistcens(d, "norm")
  fb <- fitdistcens(d, "logis")
  par(mfrow = c(1,2))
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  
  # d4 = case with also right censored data
  # with differences between the algorithms by the way they polish 
  # the ECDF function, by putting some masses to zero.
  require(actuar)
  data(fluazinam)
  d4 <- fluazinam
  d <- d4
  compare.plotdistcens(d)
  compare.npmle(d)
  fa <- fitdistcens(d, "lnorm")
  fb <- fitdistcens(d, "llogis")
  par(mfrow = c(1,2))
  cdfcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  cdfcompcens(list(fa, fb), NPMLE.method = "Wang")
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  
  
  # d5 a random example with exact values
  set.seed(123)
  r <- rnorm(10)
  d5 <- data.frame(left = r, right = r)
  d <- d5
  compare.plotdistcens(d)
  compare.npmle(d)
  fa <- fitdistcens(d, "norm")
  fb <- fitdistcens(d, "logis")
  par(mfrow = c(1,2))
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  
  # d7 = bigger dataset with also rigth censored data 
  data(salinity) 
  d7 <- log10(salinity)
  d <- d7
  compare.plotdistcens(d)
  compare.npmle(d)
  fa <- fitdistcens(d, "logis")
  fb <- fitdistcens(d, "norm")
  par(mfrow = c(1,2))
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  
  
  # d8 = an random example with all types of data (a small one)
  # set.seed(1234) # check OK
  # set.seed(1231) # check OK
  set.seed(1232)
  ns <- 25
  r <- rnorm(ns)
  d8 <- data.frame(left = r, right = r)
  delta <- rlnorm(ns)
  icensored <- rbinom(ns, size = 1, prob = 0.2) 
  Lcensored <- rbinom(ns, size = 1, prob = 0.2*(1 - icensored))
  Rcensored <- rbinom(ns, size = 1, prob = 0.3*(1 - icensored)*(1 - Lcensored))
  # icensored +  Lcensored + Rcensored
  d8$left <- d8$left * (1 - Lcensored) + (-1000) * Lcensored
  d8$right <- d8$right * (1 - Rcensored) + (1000) * Rcensored
  d8$right <- d8$right + delta * icensored
  d8$right[d8$right == 1000] <- NA
  d8$left[d8$left == -1000] <- NA
  d8
  d <- d8
  compare.plotdistcens(d)
  compare.npmle(d)
  fa <- fitdistcens(d, "logis")
  fb <- fitdistcens(d, "norm")
  par(mfrow = c(1,2))
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  ppcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Turnbull.intervals")
  qqcompcens(list(fa, fb), ynoise = FALSE, NPMLE.method = "Wang")
  
  # d8 = an random example with all types of data (a big one)
  # set.seed(1234) # check OK
  # set.seed(1231) # check OK
  set.seed(1232)
  ns <- 500
  # ns <- 5000 
  r <- rnorm(ns)
  d8 <- data.frame(left = r, right = r)
  delta <- rlnorm(ns)
  icensored <- rbinom(ns, size = 1, prob = 0.2) 
  Lcensored <- rbinom(ns, size = 1, prob = 0.2*(1 - icensored))
  Rcensored <- rbinom(ns, size = 1, prob = 0.3*(1 - icensored)*(1 - Lcensored))
  # icensored +  Lcensored + Rcensored
  d8$left <- d8$left * (1 - Lcensored) + (-1000) * Lcensored
  d8$right <- d8$right * (1 - Rcensored) + (1000) * Rcensored
  d8$right <- d8$right + delta * icensored
  d8$right[d8$right == 1000] <- NA
  d8$left[d8$left == -1000] <- NA
  d <- d8
  
  par(mfrow = c(2,2))
  system.time(plotdistcens(d, NPMLE.method = "Turnbull.middlepoints"))
  system.time(plotdistcens(d, NPMLE.method = "Turnbull.intervals"))
  system.time(plotdistcens(d, NPMLE.method = "Wang")) 
  
}

