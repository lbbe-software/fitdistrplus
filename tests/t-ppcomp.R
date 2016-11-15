library(fitdistrplus)

# ?cdfcomp


# (1) Plot various distributions fitted to serving size data
#
data(groundbeef)
serving <- groundbeef$serving
fitW <- fitdist(serving, "weibull")
fitln <- fitdist(serving, "lnorm")
fitg <- fitdist(serving, "gamma")

# sanity checks
try(ppcomp("list(fitW, fitln, fitg)"), silent = TRUE)
try(ppcomp(list(fitW, fitln, fitg, a = 1)), silent = TRUE)

# real call
ppcomp(list(fitW, fitln, fitg))

ppcomp(list(fitW, fitln, fitg), legendtext = c("Weibull", "lognormal", "gamma"),
       main="ground beef fits", xlab="Theo.",
       ylab="serving sizes (g)", xlim = c(0, 1/2))

ppcomp(list(fitW, fitln, fitg), legendtext=c("Weibull","lognormal","gamma"),
       main="ground beef fits", xlab="Theo.",
       ylab="serving sizes (g)", xlogscale=TRUE, line01=FALSE)

ppcomp(list(fitW, fitln, fitg), legendtext=c("Weibull","lognormal","gamma"),
       main="ground beef fits", xlab="Theo.",
       ylab="serving sizes (g)", ylogscale=TRUE, line01=FALSE)

ppcomp(list(fitW, fitln, fitg), legendtext=c("Weibull","lognormal","gamma"),
       main="ground beef fits", ylim=c(1e-3, 1), xlim=c(1e-3, 1),
       fitpch=c("+", "-", "."))

if (requireNamespace ("ggplot2", quietly = TRUE)) {
  ppcomp(list(fitW, fitln, fitg), plotstyle = "ggplot")
  
  ppcomp(list(fitW, fitln, fitg), legendtext=c("Weibull","lognormal","gamma"),
         main="ground beef fits", xlab="Theo.",
         ylab="serving sizes (g)", xlim = c(0, 1/2), plotstyle = "ggplot")
  
  ppcomp(list(fitW, fitln, fitg), legendtext=c("Weibull","lognormal","gamma"),
         main="ground beef fits", xlab="Theo.",
         ylab="serving sizes (g)", xlogscale=TRUE, line01=FALSE, plotstyle = "ggplot")
  
  ppcomp(list(fitW, fitln, fitg), legendtext=c("Weibull","lognormal","gamma"),
         main="ground beef fits", xlab="Theo.",
         ylab="serving sizes (g)", ylogscale=TRUE, line01=FALSE, plotstyle = "ggplot")
  
  ppcomp(list(fitW, fitln, fitg), legendtext=c("Weibull","lognormal","gamma"),
         main="ground beef fits", ylim=c(1e-3, 1), xlim=c(1e-3, 1),
         fitpch=c("+", "-", "."), plotstyle = "ggplot")
}



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

ppcomp(llfit, main="fits of a lognormal dist. using various GOF dist.")

ppcomp(llfit, xlegend="topleft",
       xlogscale=TRUE, main="fits of a lognormal dist. using various GOF dist.",
       legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"))

ppcomp(llfit, xlegend="topleft",
       xlogscale=TRUE, main="fits of a lognormal dist. using various GOF dist.",
       legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"), 
       fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"))

ppcomp(llfit, ynoise=FALSE,	xlogscale=TRUE, ylogscale=TRUE, xlim=c(1e-3,1), ylim=c(1e-3,1))

ppcomp(flnMGEKS)

if (requireNamespace ("ggplot2", quietly = TRUE)) {
  ppcomp(llfit, main="fits of a lognormal dist. using various GOF dist.", plotstyle = "ggplot")
  
  ppcomp(llfit, xlegend="topleft",
         xlogscale=TRUE, main="fits of a lognormal dist. using various GOF dist.",
         legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"), plotstyle = "ggplot")
  
  ppcomp(llfit, xlegend="topleft",
         xlogscale=TRUE, main="fits of a lognormal dist. using various GOF dist.",
         legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"), 
         fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"), plotstyle = "ggplot")
  
  ppcomp(llfit, ynoise=FALSE,	xlogscale=TRUE, ylogscale=TRUE, xlim=c(1e-3,1), ylim=c(1e-3,1), plotstyle = "ggplot")
  
  ppcomp(flnMGEKS, plotstyle = "ggplot")
}



# (3) Plot lognormal distributions fitted by 
# maximum goodness-of-fit estimation
# using various distances (data plotted in log scale)
#

x1 <- c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,13.2,8.4,6.3,8.9,5.2,10.9,14.4)
n1 <- length(x1)

dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q,a,b) exp(-exp((a-q)/b))

f1 <- mledist(x1, "norm")
f2 <- mledist(x1, "gumbel", start = list(a = 10, b = 5))
f3 <- mledist(x1, "exp")

xx1 <- pnorm(sort(x1), f1$estimate[1], f1$estimate[2])
xx2 <- pgumbel(sort(x1), f2$estimate[1], f2$estimate[2])
xx3 <- pexp(sort(x1), f3$estimate[1])
xlim <- c(xx1, xx2, xx3)
xlim <- range(xlim[which(is.finite(xlim))])

# graph 1
plot(xx1, 1:n1/n1, col = "red", xlim = xlim)
points(xx2, 1:n1/n1, col = "green")
points(xx3, 1:n1/n1, col = "blue")
legend("bottomright", pch = 1, leg = c("Normal","Gumbel","Exp"), col = c("red","green","blue"))

# graph 2
f1 <- fitdist(x1,"norm")
f2 <- fitdist(x1,"gumbel",start=list(a=10,b=5))
f3 <- fitdist(x1, "exp")
ppcomp(list(f1, f2, f3), fitcol=c("red","green","blue"), ynoise = FALSE, legendtext = c("Normal","Gumbel","Exp"))

# graph 3
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  ppcomp(list(f1, f2, f3), fitcol=c("red","green","blue"), ynoise = FALSE, legendtext = c("Normal","Gumbel","Exp"), plotstyle = "gg")
}


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

ppcomp(list(fit1, fit2, fit3), fitpch=rep(".", 3), fitcol=c("green", "red", "blue"))

if (requireNamespace ("ggplot2", quietly = TRUE)) {
  ppcomp(list(fit1, fit2, fit3), fitpch=rep(".", 3), fitcol=c("green", "red", "blue"), plotstyle = "gg")
}

# (5) large data
#

n <- 2e4
# n <- 1e2
x <- rlnorm(n)
f1 <- fitdist(x, "lnorm")
f2 <- fitdist(x, "exp")

ppcomp(list(f1, f2), fitpch=2)

if (requireNamespace ("ggplot2", quietly = TRUE)) {
  ppcomp(list(f1, f2), fitpch=2, plotstyle = "ggplot")
}



