library(fitdistrplus)
set.seed(123)
n <- 20

# (1) test qqcomp/qqcompcens on a good example
#


x <- rlnorm(n, 0, 1)
dx <- data.frame(left=x, right=x)
dx$right[1:(n/2)*2] <- NA
dx$left[2:(n/4)*4-1] <- NA

f1 <- fitdist(x, "lnorm")
f1c <- fitdistcens(dx, "lnorm")

f3 <- fitdist(x, "lnorm", fix.arg=list(sdlog=1))
f3c <- fitdistcens(dx, "lnorm", fix.arg=list(sdlog=1))



par(mfrow=1:2)
qqcomp(f1)
qqcompcens(f1c)
par(mfrow = c(1,1))
if(requireNamespace ("ggplot2", quietly = TRUE)) {
  qqcomp(f1, plotstyle = "ggplot")
  qqcompcens(f1c, plotstyle = "ggplot")
}

#test log-scale
par(mfrow=1:2, mar=c(4,4,2,1))
qqcomp(f1, xlogscale = TRUE, ylogscale = TRUE)
qqcompcens(f1c, xlogscale = TRUE, ylogscale = TRUE)
par(mfrow = c(1,1))
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  qqcomp(f1, xlogscale = TRUE, ylogscale = TRUE, plotstyle = "ggplot")
  qqcompcens(f1c, xlogscale = TRUE, ylogscale = TRUE, plotstyle = "ggplot")
}

# (2) test qqcomp/qqcompcens on a weird example
#


f2 <- fitdist(x, "unif")
f2c <- fitdistcens(dx, "unif")

par(mfrow=1:2, mar=c(4,4,2,1))
qqcomp(list(f1, f2, f3))
qqcompcens(list(f1c, f2c, f3c))
par(mfrow = c(1,1))
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  qqcomp(list(f1, f2, f3), plotstyle = "ggplot")
  qqcompcens(list(f1c, f2c, f3c), plotstyle = "ggplot")
}

#test log-scale
par(mfrow=1:2, mar=c(4,4,2,1))
qqcomp(list(f1, f2, f3), xlogscale = TRUE, ylogscale = TRUE)
qqcompcens(list(f1c, f2c, f3c), xlogscale = TRUE, ylogscale = TRUE)
par(mfrow = c(1,1))
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  qqcomp(list(f1, f2, f3), xlogscale = TRUE, ylogscale = TRUE, plotstyle = "ggplot")
  qqcompcens(list(f1c, f2c, f3c), xlogscale = TRUE, ylogscale = TRUE, plotstyle = "ggplot")
}

#test y noise
par(mfrow=1:2, mar=c(4,4,2,1))
qqcomp(list(f1, f2, f3))
qqcomp(list(f1, f2, f3), ynoise=FALSE)
par(mfrow = c(1,1))
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  qqcomp(list(f1, f2, f3), ynoise=FALSE, plotstyle = "ggplot")
}

par(mfrow=1:2, mar=c(4,4,2,1))
qqcompcens(list(f1c, f2c, f3c))
qqcompcens(list(f1c, f2c, f3c), ynoise=FALSE)
par(mfrow = c(1,1))
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  qqcompcens(list(f1c, f2c, f3c), plotstyle = "ggplot")
}

#test log-scale y-noise
par(mfrow=1:2, mar=c(4,4,2,1))
qqcomp(list(f1, f2, f3), xlogscale = TRUE, ylogscale = TRUE)
qqcomp(list(f1, f2, f3), xlogscale = TRUE, ylogscale = TRUE, ynoise=FALSE)
par(mfrow = c(1,1))
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  qqcomp(list(f1, f2, f3), xlogscale = TRUE, ylogscale = TRUE, ynoise=FALSE, plotstyle = "ggplot")
}

par(mfrow=1:2, mar=c(4,4,2,1))
qqcompcens(list(f1c, f2c), xlogscale = TRUE, ylogscale = TRUE)
qqcompcens(list(f1c, f2c), xlogscale = TRUE, ylogscale = TRUE, ynoise=FALSE)
par(mfrow = c(1,1))
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  qqcompcens(list(f1c, f2c), xlogscale = TRUE, ylogscale = TRUE, plotstyle = "ggplot")
}
