library(fitdistrplus)
set.seed(123)
n <- 20

# (1) test ppcomp/ppcompcens on a good example
#


x <- rlnorm(n, 0, 1)
dx <- data.frame(left=x, right=x)
dx$right[1:(n/2)*2] <- NA
dx$left[2:(n/4)*4-1] <- NA

f1 <- fitdist(x, "lnorm")
f1c <- fitdistcens(dx, "lnorm")

f3 <- fitdist(x, "lnorm", fix.arg=list(sdlog=1))
f3c <- fitdistcens(dx, "lnorm", fix.arg=list(sdlog=1))



par(mfrow=1:2, mar=c(4,4,2,1))
ppcomp(f1)
ppcompcens(f1c)

#test log-scale
par(mfrow=1:2, mar=c(4,4,2,1))
ppcomp(f1, xlogscale = TRUE, ylogscale = TRUE)
ppcompcens(f1c, xlogscale = TRUE, ylogscale = TRUE, xlim=c(.1, 1), ylim=c(.1, 1))

# with ggplot
ppcomp(f1, xlogscale = TRUE, ylogscale = TRUE, plotstyle = "ggplot")


# (2) test ppcomp/ppcompcens on a weird example
#


f2 <- fitdist(x, "unif")
f2c <- fitdistcens(dx, "unif")

par(mfrow=1:2, mar=c(4,4,2,1))
ppcomp(list(f1, f2, f3))
ppcompcens(list(f1c, f2c, f3c))

# with ggplot
ppcomp(list(f1, f2, f3), plotstyle = "ggplot")
#PB no noise ? Or the same for each fitted distribution (no interest)

#test log-scale
par(mfrow=1:2, mar=c(4,4,2,1))
ppcomp(list(f1, f2, f3), xlogscale = TRUE, ylogscale = TRUE, xlim=c(.1, 1), ylim=c(.1, 1))
ppcompcens(list(f1c, f2c, f3c), xlogscale = TRUE, ylogscale = TRUE, xlim=c(.1, 1), ylim=c(.1, 1))

# with ggplot
ppcomp(list(f1, f2, f3), xlogscale = TRUE, ylogscale = TRUE, 
       xlim=c(.1, 1), ylim=c(.1, 1), plotstyle = "ggplot")
#PB warning !!!!!!!!!!!!!!!!!!!!!!

#test y noise
par(mfrow=1:2, mar=c(4,4,2,1))
ppcomp(list(f1, f2, f3))
ppcomp(list(f1, f2, f3), ynoise=FALSE)

# with ggplot
ppcomp(list(f1, f2, f3), plotstyle = "ggplot")
ppcomp(list(f1, f2, f3), ynoise=FALSE, plotstyle = "ggplot")
#PB : the same noise in put for each fitted dist. -> not the purpose

par(mfrow=1:2, mar=c(4,4,2,1))
ppcompcens(list(f1c, f2c, f3c))
ppcompcens(list(f1c, f2c, f3c), ynoise=FALSE)


#test log-scale y-noise
par(mfrow=1:2, mar=c(4,4,2,1))
ppcomp(list(f1, f2, f3), xlogscale = TRUE, ylogscale = TRUE, xlim=c(.1, 1), ylim=c(.1, 1))
ppcomp(list(f1, f2, f3), xlogscale = TRUE, ylogscale = TRUE, xlim=c(.1, 1), ylim=c(.1, 1), ynoise=FALSE)

# with ggplot
ppcomp(list(f1, f2, f3), xlogscale = TRUE, ylogscale = TRUE, 
       xlim=c(.1, 1), ylim=c(.1, 1), plotstyle = "ggplot")
ppcomp(list(f1, f2, f3), xlogscale = TRUE, ylogscale = TRUE, 
       xlim=c(.1, 1), ylim=c(.1, 1), ynoise=FALSE, plotstyle = "ggplot")
#PB both pb with warning and wrong noise

par(mfrow=1:2, mar=c(4,4,2,1))
ppcompcens(list(f1c, f2c), xlogscale = TRUE, ylogscale = TRUE, xlim=c(.1, 1), ylim=c(.1, 1))
ppcompcens(list(f1c, f2c), xlogscale = TRUE, ylogscale = TRUE, xlim=c(.1, 1), ylim=c(.1, 1), ynoise=FALSE)
