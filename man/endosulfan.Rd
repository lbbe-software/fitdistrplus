\name{endosulfan}
\alias{endosulfan}
\docType{data}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Species-Sensitivity Distribution (SSD) for endosulfan }

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
(arthropods, fish or nonarthropod invertebrates).

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

# (2) plot and description of data for non Australian fish in logarithm
#
log10ATV <-log10(subset(endosulfan,(Australian == "no") & (group == "Fish"))$ATV)
plotdist(log10ATV)
descdist(log10ATV,boot=1000)

# (3) fit of a lognormal and a logistic distribution to data
# (classical distributions used for SSD)
# and visual comparison of the fits using a cdf plot
#
fln <- fitdist(log10ATV,"norm")
summary(fln)

fll <- fitdist(log10ATV,"logis")
summary(fll)

cdfcomp(list(fln,fll),legendtext=c("normal","logistic"),
xlab="log10ATV")

# (4) estimation of the 5 percent quantile value of 
# each fitted distribution (5 percent hazardous concentration  : HC5)
# with its two-sided 95 percent confidence interval calculated by 
# non parametric bootstrap
#
# in log10(ATV)
HC5ln <- quantiledist(fln,probs = 0.05,bootstrap=TRUE)
HC5ll <- quantiledist(fll,probs = 0.05,bootstrap=TRUE)
# in ATV
10^(HC5ln$quantiles)
10^(HC5ln$quantCI)
10^(HC5ll$quantiles)
10^(HC5ll$quantCI)

# (5) estimation of the 5 percent quantile value of 
# the fitted lognormal distribution (5 percent hazardous concentration  : HC5)
# with its one-sided 95 percent confidence interval (type "greater")
# calculated by non parametric bootstrap
#
# in log10(ATV)
HC5ln <- quantiledist(fln,probs = 0.05,bootstrap=TRUE,CI.type="greater")
# in ATV
10^(HC5ln$quantiles)
10^(HC5ln$quantCI)

# (6) fit of a logistic distribution 
# by minimizing the modified Anderson-Darling AD2L distance
# cf. ?mgedist for definition of this distance
#

fllAD2L <- fitdist(log10ATV,"logis",method="mge",gof="AD2L")
summary(fllAD2L)
plot(fllAD2L)

}

\keyword{ datasets }% at least one, from doc/KEYWORDS
