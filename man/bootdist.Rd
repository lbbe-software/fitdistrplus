\name{bootdist}
\alias{bootdist}
\alias{plot.bootdist}
\alias{print.bootdist}
\alias{summary.bootdist}

\title{ Bootstrap simulation of uncertainty for non-censored data}

\description{
 Uses parametric or nonparametric bootstrap resampling in order to simulate uncertainty 
 in the parameters of the distribution fitted to non-censored data.  
}

\usage{
bootdist(f, bootmethod = "param", niter = 1001, silent = TRUE, 
      parallel = c("no", "snow", "multicore"), ncpus)
\method{print}{bootdist}(x, \dots)
\method{plot}{bootdist}(x, main = "Bootstrapped values of parameters", enhance = FALSE, 
    trueval = NULL, rampcol = NULL, nbgrid = 100, nbcol = 100, \dots)
\method{summary}{bootdist}(object, \dots)

}

\arguments{
\item{f}{ An object of class \code{"fitdist"}, output of the \code{\link{fitdist}} function.}
\item{bootmethod}{ A character string coding for the type of resampling : \code{"param"} 
  for a parametric resampling and \code{"nonparam"} for a nonparametric resampling of data.}
\item{niter}{ The number of samples drawn by bootstrap.  }
\item{silent}{A logical to remove or show warnings and errors when bootstraping.}
\item{parallel}{The type of parallel operation to be used, \code{"snow"} or \code{"multicore"}
  (the second one not being available on Windows),
  or \code{"no"} if no parallel operation.}
\item{ncpus}{Number of processes to be used in parallel operation : 
  typically one would fix it to the number of available CPUs.}

\item{x}{ An object of class \code{"bootdist"}.  }
\item{object}{ An object of class \code{"bootdist"}.  }
\item{main}{an overall title for the plot: see \code{\link{title}},
default to \code{"Bootstrapped values of parameters"}.}
\item{enhance}{a logical to get an enhanced plot.}
\item{trueval}{when relevant, a numeric vector with the true value of 
parameters (for backfitting purposes).}
\item{rampcol}{colors to interpolate; must be a valid argument to 
\code{\link[grDevices]{colorRampPalette}()}.}
\item{nbgrid}{Number of grid points in each direction. Can be scalar or a length-2 integer vector.}
\item{nbcol}{an integer argument, the required number of colors}
\item{\dots}{ Further arguments to be passed to generic methods }

}

\details{
    Samples are drawn by parametric bootstrap (resampling from the distribution fitted by 
    \code{\link{fitdist}}) or nonparametric bootstrap (resampling with replacement from the 
    data set). On each bootstrap sample the function 
    \code{\link{mledist}}  (or \code{\link{mmedist}}, \code{\link{qmedist}}, \code{\link{mgedist}} 
    according to the component \code{f$method} of the object of class \code{"fitdist"}) is 
    used to estimate bootstrapped values of parameters. When that function fails
    to converge, \code{NA} values are returned. Medians and 2.5 and 97.5 percentiles are 
    computed by removing \code{NA} values. 
    The medians and the 95 percent confidence intervals of parameters (2.5 and 97.5 
    percentiles) are printed in the summary.
    If inferior to the whole number of iterations, the number of iterations for which 
    the function converges is also printed in the summary.
    
    By default (when \code{enhance=FALSE}), the plot of an object of class 
    \code{"bootdist"} consists in a scatterplot or a matrix 
    of scatterplots of the bootstrapped values of parameters.
    It uses the function \code{\link{stripchart}} when the fitted distribution
    is characterized by only one parameter, the function \code{\link{plot}} when there
    are two paramters and the function \code{\link{pairs}} in other cases.
    In these last cases, it provides a representation of the joint uncertainty distribution 
    of the fitted parameters.
    
    When \code{enhance=TRUE}, a personalized plot version of \code{\link{pairs}} is used where
    upper graphs are scatterplots and lower graphs are heatmap image using \code{\link{image}}
    based on a kernel based estimator for the 2D density function (using \code{kde2d} from
    MASS package).
    Arguments \code{rampcol}, \code{nbgrid}, \code{nbcol} can be used to customize the plots.
    Defautls values are \code{rampcol=c("green", "yellow", "orange", "red")}, \code{nbcol=100} 
    (see \code{\link[grDevices]{colorRampPalette}()}), \code{nbgrid=100} (see \code{kde2d}).
    In addition, when fitting parameters on simulated datasets for backtesting purposes, an
    additional argument \code{trueval} can be used to plot a cross at the true value.
    
    It is possible to accelerate the bootstrap using parallelization. We recommend you to
    use \code{parallel = "multicore"}, or \code{parallel = "snow"} if you work on Windows,
    and to fix \code{ncpus} to the number of available processors.
}

\value{ 
    \code{bootdist} returns an object of class \code{"bootdist"}, a list with 6 components,
    \item{estim}{ a data frame containing the bootstrapped values of parameters.}
    \item{converg}{ a vector containing the codes for convergence obtained if an iterative method is
        used to estimate parameters on 
        each bootstraped data set (and 0 if a closed formula is used).}
    \item{method}{ A character string coding for the type of resampling : \code{"param"} for a parametric
        resampling and \code{"nonparam"} for a nonparametric resampling. }
    \item{nbboot}{ The number of samples drawn by bootstrap.}   
    \item{CI}{ bootstrap medians and 95 percent confidence percentile intervals of parameters. }
    \item{fitpart}{ The object of class \code{"fitdist"} on which the bootstrap procedure
        was applied.}
        
Generic functions:
\describe{       
    \item{\code{print}}{
    The print of a \code{"bootdist"} object shows the bootstrap parameter estimates. If inferior to the whole number of bootstrap iterations, 
    the number of iterations 
    for which the estimation converges is also printed.
    }
    \item{\code{summary}}{
    The summary provides the median and 2.5 and 97.5 percentiles of each parameter. If inferior to the whole number of bootstrap iterations, 
    the number of iterations 
    for which the estimation converges is also printed in the summary.
    }
    \item{\code{plot}}{
    The plot shows the bootstrap estimates with \code{\link{stripchart}} function
    for univariate parameters and \code{\link{plot}} function for multivariate parameters.
    }           
    
}
        
}
        

\seealso{ 
    See \code{\link{fitdistrplus}} for an overview of the package.
    \code{\link{fitdist}}, \code{\link{mledist}}, \code{\link{qmedist}}, \code{\link{mmedist}},
    \code{\link{mgedist}},
    \code{\link{quantile.bootdist}} for another generic function to calculate 
    quantiles from the fitted distribution and its bootstrap results
    and \code{\link{CIcdfplot}} for adding confidence intervals on quantiles
    to a CDF plot of the fitted distribution.

}

\references{ 
Cullen AC and Frey HC (1999), \emph{Probabilistic techniques in exposure assessment}.
Plenum Press, USA, pp. 181-241.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.
}


\author{ 
Marie-Laure Delignette-Muller and Christophe Dutang.
}


\examples{
# We choose a low number of bootstrap replicates in order to satisfy CRAN running times
# constraint.
# For practical applications, we recommend to use at least niter=501 or niter=1001.


# (1) Fit of a gamma distribution to serving size data
# using default method (maximum likelihood estimation)
# followed by parametric bootstrap
#
data(groundbeef)
x1 <- groundbeef$serving
f1 <- fitdist(x1, "gamma")
b1 <- bootdist(f1, niter=51)
print(b1)
plot(b1)
plot(b1, enhance=TRUE)
summary(b1)
quantile(b1)
CIcdfplot(b1, CI.output = "quantile")

# (2) non parametric bootstrap on the same fit
#
b1b <- bootdist(f1, bootmethod="nonparam", niter=51) 
summary(b1b)
quantile(b1b)


# (3) Fit of a normal distribution on acute toxicity values of endosulfan in log10 for
# nonarthropod invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5 percent quantile value of 
# the fitted distribution, what is called the 5 percent hazardous concentration (HC5)
# in ecotoxicology, with its two-sided 95 percent confidence interval calculated by 
# parametric bootstrap
#
data(endosulfan)
ATV <- subset(endosulfan, group == "NonArthroInvert")$ATV
log10ATV <- log10(subset(endosulfan, group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV, "norm")
bln <- bootdist(fln, bootmethod = "param", niter=51)
quantile(bln, probs = c(0.05, 0.1, 0.2))

# (4) comparison of sequential and parallel versions of bootstrap
# to be tried with a greater number of iterations (1001 or more)
#
\dontrun{
niter <- 1001
data(groundbeef)
x1 <- groundbeef$serving
f1 <- fitdist(x1, "gamma")

# sequential version
ptm <- proc.time()
summary(bootdist(f1, niter = niter))
proc.time() - ptm

# parallel version using snow
require(parallel)
ptm <- proc.time()
summary(bootdist(f1, niter = niter, parallel = "snow", ncpus = 4))
proc.time() - ptm

# parallel version using multicore (not available on Windows)
ptm <- proc.time()
summary(bootdist(f1, niter = niter, parallel = "multicore", ncpus = 4))
proc.time() - ptm
}

}


\keyword{ distribution }