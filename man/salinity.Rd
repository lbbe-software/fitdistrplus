\name{salinity}
\alias{salinity}
\docType{data}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Species-Sensitivity Distribution (SSD) for salinity tolerance }

\description{
72-hour acute salinity tolerance (LC50 values)
of riverine macro-invertebrates. 
}

\usage{
data(salinity)
}
%- maybe also 'usage' for other objects documented here.

\format{ 
\code{salinity} is a data frame with 2 columns named left and right, describing
each observed LC50 value (in electrical condutivity, millisiemens per centimeter) as an interval. 
The left column contains either NA for
left censored observations, the left bound of the interval for interval censored
observations, or the observed value for non-censored observations. The right
column contains either NA for right censored observations, the right bound of
the interval for interval censored observations, or the observed value for noncensored
observations.

}

\source{ 
Kefford, B.J., Nugegoda, D., Metzeling, L., Fields, E. 2006. Validating species sensitivity
distributions using salinity tolerance of riverine macroinvertebrates in
the southern Murray-darling Basin (Vitoria, Australia). \emph{Canadian Journal of Fisheries
and Aquatic Science},
 \bold{63}, 1865-1877.
 }

%\references{  }

\examples{
# (1) load of data
#
data(salinity)

# (2) plot of data using Turnbull cdf plot
#
log10LC50 <-log10(salinity)
plotdistcens(log10LC50)

# (3) fit of a lognormal and a logistic distribution to data
# (classical distributions used for species sensitivity 
# distributions, SSD, in ecotoxicology))
# and visual comparison of the fits using Turnbull cdf plot 
#
fln <- fitdistcens(log10LC50,"norm")
summary(fln)

fll <- fitdistcens(log10LC50,"logis")
summary(fll)

cdfcompcens(list(fln,fll),legendtext=c("normal","logistic"),
    xlab = "log10(LC50)",xlim=c(0.5,2),lines01 = TRUE)

# (4) estimation of the 5 percent quantile value of 
# the normal fitted distribution (5 percent hazardous concentration  : HC5)
# with its two-sided 95 percent confidence interval calculated by 
# non parametric bootstrap
#
# in log10(LC50)
bln <- bootdistcens(fln,niter=501)
HC5ln <- quantile(bln,probs = 0.05)
# in LC50
10^(HC5ln$quantiles)
10^(HC5ln$quantCI)

# (5) estimation of the HC5 value
# with its one-sided 95 percent confidence interval (type "greater")
#
# in log10(LC50)
HC5lnb <- quantile(bln,probs = 0.05,CI.type="greater")

# in LC50
10^(HC5lnb$quantiles)
10^(HC5lnb$quantCI)

}

\keyword{ datasets }% at least one, from doc/KEYWORDS
