\name{fluazinam}
\alias{fluazinam}
\docType{data}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Species-Sensitivity Distribution (SSD) for Fluazinam }

\description{
48-hour acute toxicity values (EC50 values)
for exposure of macroinvertebrates and zooplancton to fluazinam. 
}

\usage{
data(fluazinam)
}
%- maybe also 'usage' for other objects documented here.

\format{ 
\code{fluazinam} is a data frame with 2 columns named left and right, describing
each observed EC50 value (in micrograms per liter) as an interval. 
The left column contains either NA for
left censored observations, the left bound of the interval for interval censored
observations, or the observed value for non-censored observations. The right
column contains either NA for right censored observations, the right bound of
the interval for interval censored observations, or the observed value for noncensored
observations.

}

\source{ 
Hose, G.C., Van den Brink, P.J. 2004. The species sensitivity distribution approach compared 
to a microcosm study: A case study with the fungicide fluazinam. \emph{Ecotoxicology
and Environmental Safety},
 \bold{73}, 109-122.
 }

%\references{  }

\examples{
# (1) load of data
#
data(fluazinam)

# (2) plot of data using Turnbull cdf plot
#
log10EC50 <-log10(fluazinam)
plotdistcens(log10EC50)

# (3) fit of a lognormal and a logistic distribution to data
# (classical distributions used for species sensitivity 
# distributions, SSD, in ecotoxicology))
# and visual comparison of the fits using Turnbull cdf plot 
#
fln <- fitdistcens(log10EC50,"norm")
summary(fln)

fll <- fitdistcens(log10EC50,"logis")
summary(fll)

cdfcompcens(list(fln,fll),legendtext=c("normal","logistic"),
xlab = "log10(EC50)")

# (4) estimation of the 5 percent quantile value of 
# the normal fitted distribution (5 percent hazardous concentration  : HC5)
# with its two-sided 95 percent confidence interval calculated by 
# non parametric bootstrap
#
# in log10(EC50)
bln <- bootdistcens(fln,niter=501)
HC5ln <- quantile(bln,probs = 0.05)
# in EC50
10^(HC5ln$quantiles)
10^(HC5ln$quantCI)

# (5) estimation of the 5 percent quantile value of 
# the fitted lognormal distribution (5 percent hazardous concentration  : HC5)
# with its one-sided 95 percent confidence interval (type "greater")
# calculated by non parametric bootstrap
#
# in log10(EC50)
HC5lnb <- quantile(bln,probs = 0.05,CI.type="greater")

# in LC50
10^(HC5lnb$quantiles)
10^(HC5lnb$quantCI)

}

\keyword{ datasets }% at least one, from doc/KEYWORDS
