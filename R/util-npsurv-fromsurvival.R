# ----------------------------------------------------------------------- #
#   Copyright (c) 2018 Marie Laure Delignette-Muller                      #
#                                                                         #
# Nonparametric maximum likelihood estimation from interval-censored data #
# ----------------------------------------------------------------------- #
# Extrapolation of bounds of the Turnbull intervals                       #
# form the outputs of survfit (pakage survival)                           #
# ----------------------------------------------------------------------- #

  
npsurv.fromsurvival <- function(censdata)
{
  survdata <- Surv(time = censdata$left, time2 = censdata$right, type="interval2")
  survfitted <- survfit(survdata ~ 1)
  
  # calculation of mass
  s <- survfitted$surv
  ns <- length(s)
  savant <- c(1, s[-ns])
  mass <- savant - s
  
  # calculation of bounds of Turnbull intervals (equivalence classes)
  middletime <- survfitted$time
  lem <- length(middletime)
  leftbounds <- numeric(length = lem)
  rightbounds <- numeric(length = lem)
  
  db <- censdata
  db$left[is.na(db$left)] <- -Inf
  db$right[is.na(db$right)] <- Inf
  bounds <- sort(unique(c(db$right, db$left)))
  leb <- length(bounds)
  finitebounds <- bounds[is.finite(bounds)]
  minbounds <- min(finitebounds)
  maxbounds <- max(finitebounds)
  uncensoredvalues <- censdata[censdata$left == censdata$right,]$left
  
  
  j <- 1
  for (i in 1:lem)
  {
    while(bounds[j] < middletime[i]) j <- j + 1
    if (isTRUE(all.equal(middletime[i], bounds[j])))
    {
      leftbounds[i] <- bounds[j]
      rightbounds[i] <- bounds[j]
    } else
    {
      leftbounds[i] <- bounds[j - 1]
      rightbounds[i] <- bounds[j]
    }
  }
  
  
  # correction for first and last boundsif needed
  if (!is.finite(bounds[1]) & isTRUE(all.equal(leftbounds[1], rightbounds[1])) &
      isTRUE(all.equal(leftbounds[1], minbounds)) &
      !is.element(leftbounds[1], uncensoredvalues))
  {
    leftbounds[1] <- -Inf
  }
  if (!is.finite(bounds[leb]) & isTRUE(all.equal(leftbounds[lem], rightbounds[lem])) &
      isTRUE(all.equal(leftbounds[lem], maxbounds))&
      !is.element(leftbounds[lem], uncensoredvalues))
  {
    rightbounds[lem] <- Inf
  }
  mass[lem] <- 1 - sum(mass[1:(lem - 1)])
  
  f <- data.frame(left = leftbounds, right = rightbounds, p = mass, middletime = middletime)
  # elimination of intervals withh negligible masses
  threshold <- 0.0001 / nrow(censdata)
  f <- f[f$p > threshold, ]
}
