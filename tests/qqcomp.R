library(fitdistrplus)

# ?cdfcomp


# (1) Plot various distributions fitted to serving size data
#
data(groundbeef)
serving <- groundbeef$serving
fitW <- fitdist(serving,"weibull")
fitln <- fitdist(serving,"lnorm")
fitg <- fitdist(serving,"gamma")

#sanity checks
try(qqcomp("list(fitW, fitln, fitg)"), silent=TRUE)
try(qqcomp(list(fitW, fitln, fitg, a=1)), silent=TRUE)

#real call
qqcomp(list(fitW, fitln, fitg))

qqcomp(list(fitW,fitln,fitg), legendtext=c("Weibull","lognormal","gamma"),
    main="ground beef fits", xlab="serving sizes (g)",
    ylab="F", xlim = c(0,250))

qqcomp(list(fitW,fitln,fitg), legendtext=c("Weibull","lognormal","gamma"),
    main="ground beef fits", xlab="serving sizes (g)",
    ylab="F", xlogscale=TRUE)

qqcomp(list(fitW,fitln,fitg), legendtext=c("Weibull","lognormal","gamma"),
    main="ground beef fits", xlab="serving sizes (g)",
    ylab="F", ylogscale=TRUE)

qqcomp(list(fitW,fitln,fitg), legendtext=c("Weibull","lognormal","gamma"),
    main="ground beef fits", xlab="serving sizes (g)", xlegend="bottomright",
    ylab="F", ylogscale=TRUE, ylim=c(1e-3, 1), xlim=c(0, 250))



# (2) Plot lognormal distributions fitted by 
# maximum goodness-of-fit estimation
# using various distances (data plotted in log scale)
#
data(endosulfan)
ATV <-subset(endosulfan, group == "NonArthroInvert")$ATV
flnMGEKS <- fitdist(ATV,"lnorm",method="mge",gof="KS")
flnMGEAD <- fitdist(ATV,"lnorm",method="mge",gof="AD")
flnMGEADL <- fitdist(ATV,"lnorm",method="mge",gof="ADL")
flnMGEAD2L <- fitdist(ATV,"lnorm",method="mge",gof="AD2L")
llfit <- list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L)

qqcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
	main="fits of a lognormal dist. using various GOF dist.")

qqcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
	xlogscale=TRUE,main="fits of a lognormal dist. using various GOF dist.",
	legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"))

qqcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
	xlogscale=TRUE,main="fits of a lognormal dist. using various GOF dist.",
	legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"), 
	fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"))

qqcomp(llfit, 
	xlogscale=TRUE, main="fits of a lognormal dist. using various GOF dist.",
	legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"), 
	fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"), 
	datapch="+", datacol="grey")
	
qqcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
	xlogscale=TRUE, xlim=c(10,100000), ylim=c(1e-2, 1-1e-2))


qqcomp(flnMGEKS, xlogscale=TRUE, xlim=c(10,100000))


# (3) Plot lognormal distributions fitted by 
# maximum goodness-of-fit estimation
# using various distances (data plotted in log scale)
#

x1 <- c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
x <- seq(0, 1.1*max(x1), length=100)

(f1 <- mledist(x1,"norm"))
dgumbel<-function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
(f2 <- mledist(x1,"gumbel",start=list(a=10,b=5)))
(f3 <- mledist(x1, "exp"))


pgumbel<-function(q,a,b) exp(-exp((a-q)/b))


plot(ecdf(x1))
lines(x, pnorm(x, f1$estimate[1], f1$estimate[2]), col="red")
lines(x, pgumbel(x, f2$estimate[1], f2$estimate[2]), col="green")
lines(x, pexp(x, f3$estimate[1]), col="blue")

legend("bottomright", lty=1, leg=c("Normal","Gumbel","Exp"), col=c("red","green","blue"))



f1 <- fitdist(x1,"norm")
f2 <- fitdist(x1,"gumbel",start=list(a=10,b=5))
f3 <- fitdist(x1, "exp")

cdfcomp(list(f1, f2, f3), xlim=range(x), fitcol=c("red","green","blue"))



