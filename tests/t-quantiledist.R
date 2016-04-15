library(fitdistrplus)
nbboot <- 101
# (1) Fit of a normal distribution on acute toxicity log-transformed values of endosulfan for
# nonarthropod invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5, 10 and 20 percent quantile values of 
# the fitted distribution, which are called the 5, 10, 20 percent hazardous concentrations (HC5, HC10, HC20)
# in ecotoxicology, followed with calculations of their confidence intervals with various definitions.
#
data(endosulfan)
ATV <- subset(endosulfan, group == "NonArthroInvert")$ATV
log10ATV <- log10(subset(endosulfan, group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV, "norm")
quantile(fln, probs = c(0.05, 0.1, 0.2))
bln <- bootdist(fln,niter=nbboot,bootmethod="param")
quantile(bln, probs = c(0.05, 0.1, 0.2))
quantile(bln, probs = c(0.05, 0.1, 0.2), CI.type = "greater")
quantile(bln, probs = c(0.05, 0.1, 0.2), CI.level = 0.9)

# (2) Fit of a distribution on acute salinity log-transformed tolerance 
# for riverine macro-invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5, 10 and 20 percent quantile values of 
# the fitted distribution, which are called the 5, 10, 20 percent hazardous concentrations (HC5, HC10, HC20)
# in ecotoxicology, followed with calculations of their confidence intervals with various definitions.
#
data(salinity)
log10LC50 <-log10(salinity)
flncens <- fitdistcens(log10LC50,"norm")
quantile(flncens, probs = c(0.05, 0.1, 0.2))
blncens <- bootdistcens(flncens,niter=nbboot)
quantile(blncens, probs = c(0.05, 0.1, 0.2))
quantile(blncens, probs = c(0.05, 0.1, 0.2), CI.type = "greater")
quantile(blncens, probs = c(0.05, 0.1, 0.2), CI.level = 0.9)


# (3) Estimation of quantiles of the fitted distribution (fln)
# and two-sided 95 percent confidence intervals for various 
# probabilities using non-parametric bootstrap with 101 iterations
#
bln.np <- bootdist(fln, bootmethod = "nonparam", niter = nbboot)
quantile(bln.np, probs = c(0.05, 0.1, 0.2))

# (4) Fit of a loglogistic distribution on the same acute toxicity values and
# estimation of the 5 percent quantile (HC5) of the fitted distribution 
# and associated two-sided 95 percent confidence interval 
#
fll <- fitdist(log10ATV, "logis")
bll <- bootdist(fll, bootmethod = "param", niter = nbboot)
# in log10(ATV)
HC5ll <- quantile(bll, probs = 0.05)
HC5ll
# in ATV
10^(HC5ll$basequant)
10^(HC5ll$quantCI)
