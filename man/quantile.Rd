\name{quantile}
\alias{quantile}
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
\method{quantile}{bootdist}(x, probs = seq(0.1, 0.9, by=0.1),CI.type = "two.sided", 
  CI.level = 0.95, \dots)
\method{quantile}{bootdistcens}(x, probs = seq(0.1, 0.9, by=0.1),CI.type = "two.sided", 
  CI.level = 0.95, \dots)
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
    each probability specified in \code{probs}, using the estimated parameters.
    When used with an object of class \code{"bootdist"} or \code{"bootdistcens"}, percentile 
    confidence intervals and medians etimates 
    are also calculated from the bootstrap result.
    If \code{CI.type} is \code{two.sided}, 
    the \code{CI.level} two-sided confidence intervals of quantiles are calculated.
    If \code{CI.type} is \code{less} or \code{greater}, 
    the \code{CI.level} one-sided confidence intervals of quantiles are calculated.
    The print functions show the estimated quantiles with percentile confidence intervals
    and median estimates when a bootstrap resampling has been done previously,
    and the number of bootstrap iterations 
    for which the estimation converges if it is inferior to the whole number of bootstrap iterations.
}

\value{ 
    \code{quantile} returns a list with 2 components (the first two described below) when called with an object
    of class \code{"fitdist"} or \code{"fitdistcens"} and 8 components (described below)
    when called with an object  of class
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
	\item{quantmedian}{Median of bootstrap estimates (per probability).}
    \item{CI.type}{ Type of confidence interval: either \code{"two.sided"} or one-sided 
    intervals (\code{"less"} or \code{"greater"}).}
    \item{CI.level}{ The confidence level.}
    \item{nbboot}{ The number of samples drawn by bootstrap.}
    \item{nbconverg}{ The number of iterations for which the optimization algorithm converges.}
        
}
        

\seealso{ 
    See \code{\link{fitdist}}, \code{\link{bootdist}}, \code{\link{fitdistcens}}, 
    \code{\link{bootdistcens}} and \code{\link{CIcdfplot}}.
    
    Please visit the \href{https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html}{Frequently Asked Questions}.
}


\references{

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34, \doi{https://doi.org/10.18637/jss.v064.i04}.

}

\author{ 
Marie-Laure Delignette-Muller and Christophe Dutang.
}


\examples{
# (1) Fit of a normal distribution on acute toxicity log-transformed values of 
# endosulfan for nonarthropod invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5, 10 and 20 percent quantile  
# values of the fitted distribution, which are called the 5, 10, 20 percent hazardous 
# concentrations (HC5, HC10, HC20) in ecotoxicology, followed with calculations of their
# confidence intervals with various definitions, from a small number of bootstrap 
# iterations to satisfy CRAN running times constraint.
# For practical applications, we recommend to use at least niter=501 or niter=1001.
#
data(endosulfan)
ATV <- subset(endosulfan, group == "NonArthroInvert")$ATV
log10ATV <- log10(subset(endosulfan, group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV, "norm")
quantile(fln, probs = c(0.05, 0.1, 0.2))
bln <- bootdist(fln, bootmethod="param", niter=101)
quantile(bln, probs = c(0.05, 0.1, 0.2))
quantile(bln, probs = c(0.05, 0.1, 0.2), CI.type = "greater")
quantile(bln, probs = c(0.05, 0.1, 0.2), CI.level = 0.9)

# (2) Draw of 95 percent confidence intervals on quantiles of the 
# previously fitted distribution
#
cdfcomp(fln)
q1 <- quantile(bln, probs = seq(0,1,length=101))
points(q1$quantCI[1,],q1$probs,type="l")
points(q1$quantCI[2,],q1$probs,type="l")

# (2b) Draw of 95 percent confidence intervals on quantiles of the 
# previously fitted distribution
# using the NEW function CIcdfplot
#
CIcdfplot(bln, CI.output = "quantile", CI.fill = "pink")

# (3) Fit of a distribution on acute salinity log-transformed tolerance 
# for riverine macro-invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5, 10 and 20 percent quantile
# values of the fitted distribution, which are called the 5, 10, 20 percent hazardous 
# concentrations (HC5, HC10, HC20) in ecotoxicology, followed with calculations of 
# their confidence intervals with various definitions.
# from a small number of bootstrap iterations to satisfy CRAN running times constraint.
# For practical applications, we recommend to use at least niter=501 or niter=1001.
#
data(salinity)
log10LC50 <-log10(salinity)
flncens <- fitdistcens(log10LC50,"norm")
quantile(flncens, probs = c(0.05, 0.1, 0.2))
blncens <- bootdistcens(flncens, niter = 101)
quantile(blncens, probs = c(0.05, 0.1, 0.2))
quantile(blncens, probs = c(0.05, 0.1, 0.2), CI.type = "greater")
quantile(blncens, probs = c(0.05, 0.1, 0.2), CI.level = 0.9)
}

\keyword{ distribution }