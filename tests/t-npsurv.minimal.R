library(fitdistrplus)
library(survival)

# first check of alternatives for npsurv

# to delete
# ap <- cbind(L=c(0:14,1:15), 
#             R=c(1:15, rep(Inf, 15)),
#             count=c(456, 226, 152, 171, 135, 125,  83,  74,  51,  42,  
#                     43,  34,  18,   9,   6,  39,  22,  23,  24,  107, 
#                     133, 102,  68,  64,  45,  53,  33,  27,  23,  30))
# npsurv.minimal(ap, w=1, maxit=100, tol=1e-6, verb=0, pkg="stats") 
# npsurv.minimal(ap, w=1, maxit=100, tol=1e-6, verb=0, pkg="limSolve") 

#### Comparison of plotdistcens with different NPMLE methods
if(FALSE)
{
  
# d1 = trivial case with only interval censored data
d1 <- data.frame(left = c(1, 2, 3, 4, 3, 7), right = c(2, 5, 3, 7, 8, 9))
d <- d1
par(mfrow = c(2,2))
plotdistcens(d, NPMLE.method = "Turnbull.middlepoints")
plotdistcens(d, NPMLE.method = "Turnbull.intervals")
plotdistcens(d, NPMLE.method = "Wang")

# d2 = case with left and right censored data
data(smokedfish)
d2 <- smokedfish
d <- d2
par(mfrow = c(2,2))
plotdistcens(d, NPMLE.method = "Turnbull.middlepoints")
plotdistcens(d, NPMLE.method = "Turnbull.intervals")
plotdistcens(d, NPMLE.method = "Wang")

# d3 = case with also rigth censored data
d3 <- data.frame(left = c(-1.4, 1.18, -1.4, 2, -1.4, 0),
                 right = c(1, 1.18, 2, NA, 0, 2))
d <- d3
par(mfrow = c(2,2))
plotdistcens(d, NPMLE.method = "Turnbull.middlepoints")
plotdistcens(d, NPMLE.method = "Turnbull.intervals")
plotdistcens(d, NPMLE.method = "Wang")

# d4 = case with also right censored data
# with differences between the algorithms by the way they polish 
# the ECDF function, by putting some masses to zero.
data(fluazinam)
d4 <- fluazinam
d <- d4
par(mfrow = c(2,2))
plotdistcens(d, NPMLE.method = "Turnbull.middlepoints")
plotdistcens(d, NPMLE.method = "Turnbull.intervals")
plotdistcens(d, NPMLE.method = "Wang")

# d5 a random example with exact values
set.seed(123)
r <- rnorm(50)
d5 <- data.frame(left = r, right = r)
d <- d5
par(mfrow = c(2,2))
plotdistcens(d, NPMLE.method = "Turnbull.middlepoints")
plotdistcens(d, NPMLE.method = "Turnbull.intervals")
plotdistcens(d, NPMLE.method = "Wang")

# d7 = bigger dataset with also rigth censored data 
data(salinity) 
d7 <- log10(salinity)
d <- d7
par(mfrow = c(2,2))
plotdistcens(d, NPMLE.method = "Turnbull.middlepoints")
plotdistcens(d, NPMLE.method = "Turnbull.intervals")
plotdistcens(d, NPMLE.method = "Wang")


# d8 = an random example wit all types of data
# set.seed(1234) # check OK
# set.seed(1231) # check OK
set.seed(1232)
# ns <- 1000
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
d8$right[d8$right == 1000] <- Inf
d8$left[d8$left == -1000] <- -Inf
d8
d <- d8
par(mfrow = c(2,2))
system.time(plotdistcens(d, NPMLE.method = "Turnbull.middlepoints"))
system.time(plotdistcens(d, NPMLE.method = "Turnbull.intervals"))
system.time(plotdistcens(d, NPMLE.method = "Wang"))

}

#### Comparison of cdfcompcens with different NPMLE methods
if (FALSE)
{
  
  # d1 = trivial case with only interval censored data
  d1 <- data.frame(left = c(1, 2, 3, 4, 3, 7), right = c(2, 5, 3, 7, 8, 9))
  d <- d1
  fa <- fitdistcens(d, "norm")
  fb <- fitdistcens(d, "logis")
  par(mfrow = c(2,2))
  cdfcompcens(list(fa, fb), NPMLE.method = "Turnbull.middlepoints")
  cdfcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  cdfcompcens(list(fa, fb), NPMLE.method = "Wang")
  
  # d2 = case with left and right censored data
  data(smokedfish)
  d2 <- smokedfish
  d <- d2
  fa <- fitdistcens(d, "lnorm")
  fb <- fitdistcens(d, "gamma")
  par(mfrow = c(2,2))
  cdfcompcens(list(fa, fb), NPMLE.method = "Turnbull.middlepoints", xlogscale = TRUE)
  cdfcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals", xlogscale = TRUE)
  cdfcompcens(list(fa, fb), NPMLE.method = "Wang", xlogscale = TRUE)
}

#### Comparison of PPcomp with different NPMLE methods
if(FALSE)
{
  
  # d1 = trivial case with only interval censored data
  d1 <- data.frame(left = c(1, 2, 3, 4, 3, 7), right = c(2, 5, 3, 7, 8, 9))
  d <- d1
  fa <- fitdistcens(d, "norm")
  fb <- fitdistcens(d, "logis")
  par(mfrow = c(1,2))
  ppcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  ppcompcens(list(fa, fb), NPMLE.method = "Wang")
  
  # d2 = case with left and right censored data
  data(smokedfish)
  d2 <- smokedfish
  d <- d2
  fa <- fitdistcens(d, "lnorm")
  fb <- fitdistcens(d, "gamma")
  par(mfrow = c(1,2))
  ppcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  ppcompcens(list(fa, fb), NPMLE.method = "Wang")
  
  # d4 = case with also right censored data
  # with differences between the algorithms by the way they polish 
  # the ECDF function, by putting some masses to zero.
  data(fluazinam)
  d4 <- fluazinam
  d <- d4
  fa <- fitdistcens(d, "lnorm")
  fb <- fitdistcens(d, "norm")
  par(mfrow = c(1,2))
  ppcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  ppcompcens(list(fa, fb), NPMLE.method = "Wang")
  
  # d7 = bigger dataset with also rigth censored data 
  data(salinity) 
  d7 <- log10(salinity)
  d <- d7
  fa <- fitdistcens(d, "logis")
  fb <- fitdistcens(d, "norm")
  par(mfrow = c(1,2))
  ppcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  ppcompcens(list(fa, fb), NPMLE.method = "Wang")
  
}

#### Comparison of QQcomp with different NPMLE methods
if(FALSE)
{
  
  # d1 = trivial case with only interval censored data
  d1 <- data.frame(left = c(1, 2, 3, 4, 3, 7), right = c(2, 5, 3, 7, 8, 9))
  d <- d1
  fa <- fitdistcens(d, "norm")
  fb <- fitdistcens(d, "logis")
  par(mfrow = c(1,2))
  qqcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  qqcompcens(list(fa, fb), NPMLE.method = "Wang")
  
  # d2 = case with left and right censored data
  data(smokedfish)
  d2 <- smokedfish
  d <- d2
  fa <- fitdistcens(d, "lnorm")
  fb <- fitdistcens(d, "gamma")
  par(mfrow = c(1,2))
  qqcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  qqcompcens(list(fa, fb), NPMLE.method = "Wang")
  
  # d4 = case with also right censored data
  # with differences between the algorithms by the way they polish 
  # the ECDF function, by putting some masses to zero.
  data(fluazinam)
  d4 <- fluazinam
  d <- d4
  fa <- fitdistcens(d, "lnorm")
  fb <- fitdistcens(d, "norm")
  par(mfrow = c(1,2))
  qqcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  qqcompcens(list(fa, fb), NPMLE.method = "Wang")
  
  # d7 = bigger dataset with also rigth censored data 
  data(salinity) 
  d7 <- log10(salinity)
  d <- d7
  fa <- fitdistcens(d, "logis")
  fb <- fitdistcens(d, "norm")
  par(mfrow = c(1,2))
  qqcompcens(list(fa, fb), NPMLE.method = "Turnbull.intervals")
  qqcompcens(list(fa, fb), NPMLE.method = "Wang")
  
}
