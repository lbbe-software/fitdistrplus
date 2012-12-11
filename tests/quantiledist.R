library(fitdistrplus)

# (1) Fit of a lognormal distribution on acute toxicity values of endosulfan for
# nonarthropod invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5 percent quantile value of 
# the fitted distribution, what is called the 5 percent hazardous concentration (HC5)
# in ecotoxicology, with its two-sided 95 percent confidence interval calculated by 
# parametric bootstrap
#
data(endosulfan)
ATV <-subset(endosulfan,group == "NonArthroInvert")$ATV
log10ATV <-log10(subset(endosulfan,group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV, "norm")
quantile(fln, probs=c(.05, .1))
# in log10(ATV)
bln <- bootdist(fln, bootmethod = "param", niter = 101)
HC5ln <- quantile(bln, probs = 0.05)
HC5ln
# in ATV
10^(HC5ln$basequant)
10^(HC5ln$quantCI)

quantile(bln, probs = 0.05, CI.level=.85)


# (2) Estimation of quantiles of the same fitted distribution 
# and two-sided 95 percent confidence intervals for various probabilities 
#
quantile(bln, probs = c(0.05, 0.1, 0.2))


# (3) Estimation of quantiles of the same fitted distribution 
# and one-sided 95 percent confidence intervals (type "greater") for various 
# probabilities 
#
quantile(bln, probs = c(0.05, 0.1, 0.2), CI.type = "greater")

# (4) Estimation of quantiles of the fitted distribution 
# and two-sided 95 percent confidence intervals for various 
# probabilities using non-parametric bootstrap with 101 iterations
#
bln.np <- bootdist(fln, bootmethod = "nonparam", niter = 101)
quantile(bln.np, probs = c(0.05, 0.1, 0.2))

# (5) Fit of a loglogistic distribution on the same acute toxicity values and
# estimation of the 5 percent quantile (HC5) of the fitted distribution 
# and associated two-sided 95 percent confidence interval 
#
fll <- fitdist(log10ATV, "logis")
bll <- bootdist(fll, bootmethod = "param", niter = 101)
# in log10(ATV)
HC5ll <- quantile(bll, probs = 0.05)
HC5ll
# in ATV
10^(HC5ll$basequant)
10^(HC5ll$quantCI)
