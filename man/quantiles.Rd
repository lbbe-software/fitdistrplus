\name{quantiles}
\alias{quantile.fitdist}
\alias{quantile.fitdistcens}
\alias{quantile.bootdist}
\alias{quantile.bootdistcens}
\alias{print.quantile.fitdist}
\alias{print.quantile.fitdistcens}
\alias{print.quantile.bootdist}
\alias{print.quantile.bootdistcens}

\title{Quantile estimation from a fitted distribution}

\description{
 Quantile estimation from a fitted distribution, optionally with confidence intervals calculated 
 from the bootstrap result.  
}

\usage{
\method{quantile}{fitdist}(x, probs = seq(0.1, 0.9, by=0.1), \dots)
\method{quantile}{fitdistcens}(x, probs = seq(0.1, 0.9, by=0.1), \dots)
\method{quantile}{bootdist}(x, probs = seq(0.1, 0.9, by=0.1),CI.type = "two.sided", CI.level = 0.95, \dots)
\method{quantile}{bootdistcens}(x, probs = seq(0.1, 0.9, by=0.1),CI.type = "two.sided", CI.level = 0.95, \dots)
\method{print}{quantile.fitdist}(x, \dots)
\method{print}{quantile.fitdistcens}(x, \dots)
\method{print}{quantile.bootdist}(x, \dots)
\method{print}{quantile.bootdistcens}(x, \dots)

}

\arguments{
\item{x}{ An object of class \code{"fitdist"}, \code{"fitdistcens"}, \code{"bootdist"}, \code{"bootdistcens"}
    or \code{"quantile.fitdist"}, \code{"quantile.fitdistcens"}, \code{"quantile.bootdist"}, 
    \code{"quantile.bootdistcens"} for the \code{print} generic function.  }
\item{probs}{ A numeric vector of probabilities with values in [0, 1] at which quantiles must be calculated.}
\item{CI.type}{ Type of confidence intervals : either \code{"two.sided"} or one-sided 
    intervals (\code{"less"} or \code{"greater"}).}
\item{CI.level}{ The confidence level.}
\item{\dots}{Further arguments to be passed to generic functions.}
}

\details{
    Quantiles of the parametric distribution are calculated for 
    each probability specified in \code{probs}, using the estimated parameters of \code{fitpart}.
    When used with an object of class \code{"bootdist"} or \code{"bootdistcens"}, confidence intervals 
    are also calculated from the bootstrap result.
    If \code{CI.type} is \code{two.sided}, 
    the \code{CI.level} two-sided confidence intervals of quantiles are calculated.
    If \code{CI.type} is \code{less} or \code{greater}, 
    the \code{CI.level} one-sided confidence intervals of quantiles are calculated.
    The print functions show the estimated quantiles with confidence intervals if calculated,
    and the number of bootstrap iterations 
    for which the estimation converges if it is inferior to the whole number of bootstrap iterations.
}

\value{ 
    \code{quantile} returns a list with 2 components (the first two) when called with an object
    of class \code{"fitdist"} or \code{"fitdistcens"} and 8 components when called with an object  of class
    \code{"bootdist"} or \code{"bootdistcens"} :
    \item{quantiles}{ a dataframe containing the estimated quantiles for each probability value specified in 
        the argument \code{probs} (one row, and as many columns as values in \code{probs}).}
    \item{probs}{the numeric vector of probabilities at which quantiles are calculated.}
    \item{bootquant }{ A data frame containing the bootstraped values for each quantile
         (many rows, as specified in the call to \code{\link{bootdist}}  in  the argument \code{niter}, 
         and as many columns as values in \code{probs})}
    \item{quantCI}{ If \code{CI.type} is \code{two.sided}, the two bounds of the  \code{CI.level} percent 
        two.sided confidence interval for each quantile 
        (two rows and as many columns as values in \code{probs}). If \code{CI.type} is \code{less}, 
        right bound of the \code{CI.level} percent
        one.sided confidence interval for each quantile (one row).
        If \code{CI.type} is \code{greater}, left bound of the \code{CI.level} percent
        one.sided confidence interval for each quantile (one row). }
    \item{CI.type}{ Type of confidence interval: either \code{"two.sided"} or one-sided 
    intervals (\code{"less"} or \code{"greater"}).}
    \item{CI.level}{ The confidence level.}
    \item{nbboot}{ The number of samples drawn by bootstrap.}
    \item{nbconverg}{ The number of iterations fro which the optimization algorithm converges.}
        
}
        

\seealso{ 
    \code{\link{fitdist}}, \code{\link{bootdist}}, \code{\link{fitdistcens}} and \code{\link{bootdistcens}}.
}

\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} 
}


\examples{
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
bln <- bootdist(fln,niter=501,bootmethod="param")
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
blncens <- bootdistcens(flncens,niter=501)
quantile(blncens, probs = c(0.05, 0.1, 0.2))
quantile(blncens, probs = c(0.05, 0.1, 0.2), CI.type = "greater")
quantile(blncens, probs = c(0.05, 0.1, 0.2), CI.level = 0.9)


}


\keyword{ distribution }% at least one, from doc/KEYWORDS
