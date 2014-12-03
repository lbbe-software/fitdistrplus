\name{endosulfan}
\alias{endosulfan}
\docType{data}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Species Sensitivity Distribution (SSD) for endosulfan }

\description{
Summary of 48- to 96-hour acute toxicity values (LC50 and EC50 values)
for exposure of Australian an Non-Australian taxa to endosulfan. 
}

\usage{
data(endosulfan)
}
%- maybe also 'usage' for other objects documented here.

\format{ 
\code{endosulfan} is a data frame with 3 columns, named ATV for Acute Toxicity Value
(geometric mean of LC50 ou EC50 values in micrograms per liter),
Australian (coding for Australian or another origin) and group 
(arthropods, fish or non-arthropod invertebrates).

}

\source{ 
Hose, G.C., Van den Brink, P.J. 2004. Confirming the Species-Sensitivity Distribution Concept
for Endosulfan Using Laboratory, Mesocosms, and Field Data. \emph{Archives of Environmental 
Contamination and Toxicology},
 \bold{47}, 511-520.
 }

%\references{  }

\examples{
# (1) load of data
#
data(endosulfan)

# (2) plot and description of data for non Australian fish in decimal logarithm
#
log10ATV <-log10(subset(endosulfan,(Australian == "no") & (group == "Fish"))$ATV)
plotdist(log10ATV)
descdist(log10ATV,boot=1000)

# (3) fit of a normal and a logistic distribution to data in log10
# (classical distributions used for SSD)
# and visual comparison of the fits
#
fln <- fitdist(log10ATV,"norm")
summary(fln)

fll <- fitdist(log10ATV,"logis")
summary(fll)

cdfcomp(list(fln,fll),legendtext=c("normal","logistic"),
xlab="log10ATV")

denscomp(list(fln,fll),legendtext=c("normal","logistic"),
xlab="log10ATV")

qqcomp(list(fln,fll),legendtext=c("normal","logistic"))
ppcomp(list(fln,fll),legendtext=c("normal","logistic"))

gofstat(list(fln,fll), fitnames = c("lognormal", "loglogistic"))

# (4) estimation of the 5 percent quantile value of 
# logistic fitted distribution (5 percent hazardous concentration  : HC5)
# with its two-sided 95 percent confidence interval calculated by 
# parametric bootstrap 
# with a small number of iterations to satisfy CRAN running times constraint.
# For practical applications, we recommend to use at least niter=501 or niter=1001.
#
# in log10(ATV)
bll <- bootdist(fll,niter=101)
HC5ll <- quantile(bll,probs = 0.05)
# in ATV
10^(HC5ll$quantiles)
10^(HC5ll$quantCI)

# (5) estimation of the 5 percent quantile value of 
# the fitted logistic distribution (5 percent hazardous concentration  : HC5)
# with its one-sided 95 percent confidence interval (type "greater")
# calculated by 
# nonparametric bootstrap 
# with a small number of iterations to satisfy CRAN running times constraint.
# For practical applications, we recommend to use at least niter=501 or niter=1001.
# 
# in log10(ATV)
bllnonpar <- bootdist(fll,niter=101,bootmethod = "nonparam")
HC5llgreater <- quantile(bllnonpar,probs = 0.05, CI.type="greater")
# in ATV
10^(HC5llgreater$quantiles)
10^(HC5llgreater$quantCI)

# (6) fit of a logistic distribution 
# by minimizing the modified Anderson-Darling AD2L distance
# cf. ?mgedist for definition of this distance
#

fllAD2L <- fitdist(log10ATV,"logis",method="mge",gof="AD2L")
summary(fllAD2L)
plot(fllAD2L)

}

\keyword{ datasets }% at least one, from doc/KEYWORDS
