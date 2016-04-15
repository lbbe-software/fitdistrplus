library(fitdistrplus)

# ?denscomp


# (1) Plot various distributions fitted to serving size data
#
data(groundbeef)
serving <- groundbeef$serving
fitW <- fitdist(serving,"weibull")
fitln <- fitdist(serving,"lnorm")
fitg <- fitdist(serving,"gamma")

#sanity checks
try(denscomp("list(fitW, fitln, fitg)",horizontals = FALSE), silent=TRUE)
try(denscomp(list(fitW, fitln, fitg, a=1),horizontals = FALSE), silent=TRUE)

#real call
denscomp(list(fitW, fitln, fitg), probability=TRUE)
denscomp(list(fitW, fitln, fitg), probability=FALSE)

#test ylim argument
denscomp(list(fitW, fitln, fitg), probability=TRUE, ylim=c(0, .05))
denscomp(list(fitW, fitln, fitg), probability=FALSE, ylim=c(0, 100))

#test xlim, legend, main, demp    
denscomp(list(fitW, fitln, fitg), legendtext=c("Weibull","lognormal","gamma"),
    main="ground beef fits",xlab="serving sizes (g)",
    ylab="F",xlim = c(0,250), xlegend = "topright", demp=TRUE)



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

denscomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
	main="fits of a lognormal dist. using various GOF dist.")

denscomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
	main="fits of a lognormal dist. using various GOF dist.",
	legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"))

denscomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
	main="fits of a lognormal dist. using various GOF dist.",
	legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"), 
	fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"))

denscomp(llfit,
	main="fits of a lognormal dist. using various GOF dist.",
	legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"), 
	fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"), 
	datacol="grey")
	


denscomp(flnMGEKS, xlim=c(10,100000))



# (3) 
#
#

x1 <- c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
x <- seq(0, 1.1*max(x1), length=100)

f1 <- mledist(x1,"norm")
dgumbel<-function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
f2 <- mledist(x1,"gumbel",start=list(a=10,b=5))
f3 <- mledist(x1, "exp")

hist(x1, 10, prob=TRUE)
lines(x, dnorm(x, f1$estimate[1], f1$estimate[1]), col="red")
lines(x, dgumbel(x, f2$estimate[1], f2$estimate[1]), col="green")
lines(x, dexp(x, f3$estimate[1]), col="blue")

legend("topright", lty=1, leg=c("Normal","Gumbel","Exp"), col=c("red","green","blue"))

pgumbel<-function(x,a,b) exp(-exp((a-x)/b))
f1 <- fitdist(x1,"norm")
f2 <- fitdist(x1,"gumbel",start=list(a=10,b=5))
f3 <- fitdist(x1, "exp")

denscomp(list(f1, f2, f3), xlim=c(0, 30), ylim=c(0, .1))


# (4) normal mixture
#

#mixture of two normal distributions
#density
dnorm2 <- function(x, poid, m1, s1, m2, s2)
	poid*dnorm(x, m1, s1) + (1-poid)*dnorm(x, m2, s2)
#numerical approximate quantile function
qnorm2 <- function(p, poid, m1, s1, m2, s2)
{
	L2 <- function(x, prob)
		(prob - pnorm2(x, poid, m1, s1, m2, s2))^2	
	sapply(p, function(pr) optimize(L2, c(-1000, 1000), prob=pr)$minimum)
}	
#distribution function		
pnorm2 <- function(q, poid, m1, s1, m2, s2)
	poid*pnorm(q, m1, s1) + (1-poid)*pnorm(q, m2, s2)		


#basic normal distribution
set.seed(1234)
x2 <- c(rnorm(1000, 5),  rnorm(1000, 10))
#MLE fit
fit1 <- fitdist(x2, "norm2", "mle", start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2), 
	lower=c(0, 0, 0, 0, 0))
fit2 <- fitdist(x2, "norm2", "qme", probs=c(1/6, 1/4, 1/3, 1/2, 2/3), 
	start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2), 
	lower=c(0, 0, 0, 0, 0), upper=c(1/2, Inf, Inf, Inf, Inf))
fit3 <- fitdist(x2, "norm2", "mge", gof="AD", 
	start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2), 
	lower=c(0, 0, 0, 0, 0), upper=c(1/2, Inf, Inf, Inf, Inf))

denscomp(list(fit1, fit2, fit3), datapch=".")


# (5) large data
#

n <- 2e4
n <- 1e2
x <- rnorm(n)
f <- fitdist(x, "norm")

denscomp(f)
denscomp(f, demp=TRUE)



