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
fln <- fitdist(log10ATV,"norm")
# in log10(ATV)
HC5ln <- quantile(fln,probs = 0.05,bootstrap=TRUE, 
	bootstrap.arg = list(bootmethod = "param",niter = 101))
HC5ln$quantCI
# in ATV
10^(HC5ln$quantiles)
10^(HC5ln$quantCI)


# (2) Estimation of quantiles of the same fitted distribution 
# and two-sided 95 percent confidence intervals for various probabilities 
#
quantile(fln, probs = c(0.05,0.1,0.2), bootstrap = TRUE)


# (3) Estimation of quantiles of the same fitted distribution 
# and one-sided 95 percent confidence intervals (type "greater") for various 
# probabilities 
#
quantile(fln, probs = c(0.05,0.1,0.2), bootstrap = TRUE,CI.type = "greater")

# (4) Estimation of quantiles of the fitted distribution 
# and two-sided 95 percent confidence intervals for various 
# probabilities using non-parametric bootstrap with 101 iterations
#
quantile(fln,probs = c(0.05,0.1,0.2), bootstrap = TRUE,
    bootstrap.arg = list(bootmethod = "nonparam", niter = 101))

# (5) Fit of a loglogistic distribution on the same acute toxicity values and
# estimation of the 5 percent quantile (HC5) of the fitted distribution 
# and associated two-sided 95 percent confidence interval 
#
fll <- fitdist(log10ATV,"logis")
# in log10(ATV)
HC5ll <- quantile(fll,probs = 0.05,bootstrap=TRUE, 
	bootstrap.arg = list(bootmethod = "param",niter = 101))
HC5ll
# in ATV
10^(HC5ll$quantiles)
10^(HC5ll$quantCI)
