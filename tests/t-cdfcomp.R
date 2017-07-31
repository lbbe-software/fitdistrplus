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
try(cdfcomp("list(fitW, fitln, fitg)",horizontals = FALSE), silent=TRUE)
try(cdfcomp(list(fitW, fitln, fitg, a=1),horizontals = FALSE), silent=TRUE)

#real call
cdfcomp(list(fitW, fitln, fitg), horizontals = FALSE)
cdfcomp(list(fitW, fitln, fitg), horizontals = TRUE)
cdfcomp(list(fitW, fitln, fitg), horizontals = TRUE, lines01 = TRUE)
cdfcomp(list(fitW, fitln, fitg), horizontals = TRUE, verticals = TRUE, datacol = "grey")
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(fitW, fitln, fitg), horizontals = FALSE, plotstyle = "ggplot")
  cdfcomp(list(fitW, fitln, fitg), horizontals = TRUE, plotstyle = "ggplot")
  cdfcomp(list(fitW, fitln, fitg), horizontals = TRUE, lines01 = TRUE, plotstyle = "ggplot")
  cdfcomp(list(fitW, fitln, fitg), horizontals = TRUE, verticals = TRUE, datacol = "grey", plotstyle = "ggplot")
}

cdfcomp(list(fitW, fitln, fitg), legendtext = c("Weibull", "lognormal", "gamma"),
        main = "ground beef fits", xlab = "serving sizes (g)",
        ylab = "F(g)", xlim = c(0, 250), ylim = c(.5, 1))
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(fitW, fitln, fitg), legendtext = c("Weibull", "lognormal", "gamma"),
          main = "ground beef fits", xlab = "serving sizes (g)",
          ylab = "F(g)", xlim = c(0, 250), ylim = c(.5, 1), plotstyle = "ggplot")
}

cdfcomp(list(fitW, fitln, fitg), legendtext = c("Weibull", "lognormal", "gamma"),
        main = "ground beef fits", xlab = "serving sizes (g)",
        ylab = "F(g)", xlogscale = TRUE)
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(fitW, fitln, fitg), legendtext = c("Weibull", "lognormal", "gamma"),
          main = "ground beef fits", xlab = "serving sizes (g)",
          ylab = "F(g)", xlogscale = TRUE, plotstyle = "ggplot")
}

cdfcomp(list(fitW,fitln,fitg),legendtext=c("Weibull","lognormal","gamma"),
        main="ground beef fits",xlab="serving sizes (g)",
        ylab="F(g)", xlogscale=TRUE, ylogscale=TRUE, ylim=c(.005, .99))
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(fitW,fitln,fitg),legendtext=c("Weibull","lognormal","gamma"),
          main="ground beef fits",xlab="serving sizes (g)",
          ylab="F(g)", xlogscale=TRUE, ylogscale=TRUE, ylim=c(.005, .99), plotstyle = "ggplot")
}

cdfcomp(list(fitW,fitln,fitg), legendtext=c("Weibull","lognormal","gamma"),
        main="ground beef fits",xlab="serving sizes (g)",
        ylab="F(g)",xlim = c(0,250), xlegend = "topleft")
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(fitW,fitln,fitg), legendtext=c("Weibull","lognormal","gamma"),
          main="ground beef fits",xlab="serving sizes (g)",
          ylab="F(g)",xlim = c(0,250), xlegend = "topleft", plotstyle = "ggplot")
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

cdfcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L), xlogscale=TRUE, 
        main="fits of a lognormal dist. using various GOF dist.")
cdfcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
        xlogscale=TRUE,main="fits of a lognormal dist. using various GOF dist.",
        legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"))
cdfcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
        xlogscale=TRUE,main="fits of a lognormal dist. using various GOF dist.",
        legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"), 
        fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"), 
        horizontals=FALSE, datapch="+")
cdfcomp(llfit, xlogscale=TRUE, main="fits of a lognormal dist. using various GOF dist.",
        legendtext=paste("MGE", c("KS","AD","ADL","AD2L")), fitcol="grey35", 	
        fitlty="dotted", horizontals=FALSE, datapch=21, datacol="grey30")
cdfcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
        xlogscale=TRUE, verticals=TRUE, xlim=c(10,100000), datapch=21)
cdfcomp(flnMGEKS, xlogscale=TRUE, verticals=TRUE, xlim=c(10,100000))

if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L), xlogscale=TRUE,
          main = "fits of a lognormal dist. using various GOF dist.", plotstyle = "ggplot")
  cdfcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
          xlogscale=TRUE,main="fits of a lognormal dist. using various GOF dist.",
          legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"), plotstyle = "ggplot")
  cdfcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
          xlogscale=TRUE,main="fits of a lognormal dist. using various GOF dist.",
          legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"), 
          fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"), 
          horizontals=FALSE, datapch="+", plotstyle = "ggplot")
  cdfcomp(llfit, xlogscale=TRUE, main="fits of a lognormal dist. using various GOF dist.",
          legendtext=paste("MGE", c("KS","AD","ADL","AD2L")), fitcol="grey35", 	
          fitlty="dotted", horizontals=FALSE, datapch=21, datacol="grey30", plotstyle = "ggplot")
  cdfcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
          xlogscale=TRUE, verticals=TRUE, xlim=c(10,100000), datapch=21, plotstyle = "ggplot")
  cdfcomp(flnMGEKS, xlogscale=TRUE, verticals=TRUE, xlim=c(10,100000), plotstyle = "ggplot")
}




# (3) Plot normal and logistic distributions fitted by 
# maximum likelihood estimation
# using various plotting positions in cdf plots
#
data(endosulfan)
log10ATV <-log10(subset(endosulfan, group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV, "norm")
fll <- fitdist(log10ATV, "logis")

# default plot using Hazen plotting position: (1:n - 0.5)/n
cdfcomp(list(fln, fll), legendtext = c("normal", "logistic"), xlab = "log10ATV")
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(fln, fll), legendtext = c("normal", "logistic"),xlab = "log10ATV", plotstyle = "ggplot")
}

# plot using mean plotting position (named also Gumbel plotting position)
# (1:n)/(n + 1)
cdfcomp(list(fln,fll),legendtext=c("normal","logistic"),xlab="log10ATV",
        use.ppoints = TRUE, a.ppoints = 0)
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(fln,fll),legendtext=c("normal","logistic"),xlab="log10ATV",
          use.ppoints = TRUE, a.ppoints = 0, plotstyle = "ggplot")
}

# plot using basic plotting position: (1:n)/n
cdfcomp(list(fln,fll),legendtext=c("normal","logistic"),xlab="log10ATV",
        use.ppoints = FALSE)
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(fln,fll),legendtext=c("normal","logistic"),xlab="log10ATV",
          use.ppoints = FALSE, plotstyle = "ggplot")
}




# (4) Plot lognormal distributions fitted by 
# maximum goodness-of-fit estimation
# using various distances (data plotted in log scale)
#
x1 <- c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,13.2,8.4,6.3,8.9,5.2,10.9,14.4)
x <- seq(0, 1.1*max(x1), length=100)

dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q,a,b) exp(-exp((a-q)/b))

f1 <- mledist(x1, "norm")
f2 <- mledist(x1, "gumbel", start = list(a = 10, b = 5))
f3 <- mledist(x1, "exp")

# graph 1
plot(ecdf(x1))
lines(x, pnorm(x, f1$estimate[1], f1$estimate[2]), col = "red")
lines(x, pgumbel(x, f2$estimate[1], f2$estimate[2]), col = "green")
lines(x, pexp(x, f3$estimate[1]), col = "blue")
legend("bottomright", lty = 1, leg = c("Normal", "Gumbel", "Exp"), col = c("red", "green", "blue"))

# graph 2
f1 <- fitdist(x1, "norm")
f2 <- fitdist(x1, "gumbel", start = list(a = 10, b = 5))
f3 <- fitdist(x1, "exp")
cdfcomp(list(f1, f2, f3), xlim=range(x), fitcol = c("red", "green", "blue"), fitlty = 1, legendtext = c("Normal", "Gumbel", "Exp"))

# graph 3
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(f1, f2, f3), xlim=range(x), fitcol=c("red","green","blue"), fitlty = 1, legendtext = c("Normal", "Gumbel", "Exp"), plotstyle = "ggplot")
}




# (5) normal mixture
#

# mixture of two normal distributions
# density
dnorm2 <- function(x, poid, m1, s1, m2, s2)
  poid*dnorm(x, m1, s1) + (1-poid)*dnorm(x, m2, s2)
# numerically approximate quantile function
qnorm2 <- function(p, poid, m1, s1, m2, s2)
{
  L2 <- function(x, prob)
    (prob - pnorm2(x, poid, m1, s1, m2, s2))^2	
  sapply(p, function(pr) optimize(L2, c(-1000, 1000), prob=pr)$minimum)
}	
# distribution function		
pnorm2 <- function(q, poid, m1, s1, m2, s2)
  poid*pnorm(q, m1, s1) + (1-poid)*pnorm(q, m2, s2)		


# basic normal distribution
set.seed(1234)
x2 <- c(rnorm(1000, 5),  rnorm(1000, 10))
# MLE fit
fit1 <- fitdist(x2, "norm2", "mle", start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2), 
                lower=c(0, 0, 0, 0, 0))
fit2 <- fitdist(x2, "norm2", "qme", probs=c(1/6, 1/4, 1/3, 1/2, 2/3), 
                start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2), 
                lower=c(0, 0, 0, 0, 0), upper=c(1/2, Inf, Inf, Inf, Inf))
fit3 <- fitdist(x2, "norm2", "mge", gof="AD", 
                start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2), 
                lower=c(0, 0, 0, 0, 0), upper=c(1/2, Inf, Inf, Inf, Inf))

cdfcomp(list(fit1, fit2, fit3), datapch=".")
cdfcomp(list(fit1, fit2, fit3), datapch=".", xlim=c(6, 8), ylim=c(.4, .55))

if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(fit1, fit2, fit3), datapch=".", plotstyle = "ggplot")
  cdfcomp(list(fit1, fit2, fit3), datapch=".", xlim=c(6, 8), ylim=c(.4, .55), plotstyle = "ggplot")
}




# (6) discrete example
#
set.seed(1234)
x3 <- rpois(20, 10)

fit1 <- fitdist(x3, "pois", "mle")
fit2 <- fitdist(x3, "nbinom", "qme", probs=c(1/3, 2/3))

cdfcomp(list(fit1, fit2), datapch=21, horizontals=FALSE)
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(list(fit1, fit2), datapch=21, horizontals=FALSE, plotstyle = "ggplot")
}


# (7) large dataset
#
n <- 2e4
n <- 1e2
f1 <- fitdist(rlnorm(n), "lnorm")

cdfcomp(f1, do.points=TRUE)
cdfcomp(f1, do.points=FALSE)
cdfcomp(f1, horizontals = FALSE, verticals = FALSE, do.points = FALSE)
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  cdfcomp(f1, do.points = TRUE, plotstyle = "ggplot")
  cdfcomp(f1, do.points = FALSE, plotstyle = "ggplot")
  cdfcomp(f1, horizontals = FALSE, verticals = FALSE, do.points = FALSE, plotstyle = "ggplot")
}

# (8) argument add (must give the same plot (except colors) as ex. 6)
#
set.seed(1234)
x3 <- rpois(20, 10)

fit1 <- fitdist(x3, "pois", "mle")
cdfcomp(fit1, fitcol = "red", horizontals=FALSE, addlegend = FALSE)
fit2 <- fitdist(x3, "nbinom", "qme", probs=c(1/3, 2/3))
cdfcomp(fit2, fitcol = "blue", horizontals=FALSE, addlegend = FALSE, add = TRUE)
cdfcomp(list(fit1, fit2), horizontals=FALSE, addlegend = FALSE, 
        fitcol=c("red", "blue"))


if (requireNamespace ("ggplot2", quietly = TRUE)) {
  # the argument add is not available when plotstyle = "ggplot"
  cdfcomp(list(fit1, fit2), fitcol = c("red", "blue"), fitlty = 1, horizontals = FALSE, addlegend = FALSE, plotstyle = "ggplot")
}


# (9) test legend labels
#
serving <- groundbeef$serving
fitW <- fitdist(serving,"weibull")
fitW2 <- fitdist(serving,"weibull", method="qme", probs=c(1/3,2/3))
fitW3 <- fitdist(serving,"weibull", method="qme", probs=c(1/2,2/3))
fitln <- fitdist(serving,"lnorm")
fitg <- fitdist(serving,"gamma")

cdfcomp(list(fitW, fitln, fitg)) #distrib
cdfcomp(list(fitW, fitW2, fitln, fitg)) #distrib+method
cdfcomp(list(fitW, fitW2, fitW3, fitln, fitg)) #distrib+method+num
if (requireNamespace ("ggplot2", quietly = TRUE))
  cdfcomp(list(fitW, fitW2, fitW3, fitln, fitg), plotstyle = "ggplot") #distrib+method+num
