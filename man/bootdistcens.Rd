\name{bootdistcens}
\alias{bootdistcens}
\alias{plot.bootdistcens}
\alias{print.bootdistcens}
\alias{summary.bootdistcens}

\title{ Bootstrap simulation of uncertainty for censored data}

\description{
 Uses nonparametric bootstrap resampling in order to simulate uncertainty 
 in the parameters of the distribution fitted to censored data.  
}

\usage{
bootdistcens(f, niter = 1001, silent = TRUE, 
      parallel = c("no", "snow", "multicore"), ncpus)
\method{print}{bootdistcens}(x, \dots)
\method{plot}{bootdistcens}(x, \dots)
\method{summary}{bootdistcens}(object, \dots)

}

\arguments{
\item{f}{ An object of class \code{"fitdistcens"}, output of the \code{\link{fitdistcens}} function.}
\item{niter}{ The number of samples drawn by bootstrap.}
\item{silent}{A logical to remove or show warnings and errors when bootstraping.}
\item{parallel}{The type of parallel operation to be used, \code{"snow"} or \code{"multicore"}
  (the second one not being available on Windows),
  or \code{"no"} if no parallel operation.}
\item{ncpus}{Number of processes to be used in parallel operation : 
  typically one would fix it to the number of available CPUs.}

\item{x}{ An object of class \code{"bootdistcens"}.}
\item{object}{ An object of class \code{"bootdistcens"}.}
\item{\dots}{ Further arguments to be passed to generic methods.}
}

\details{
    Samples are drawn by 
    nonparametric bootstrap (resampling with replacement from the data set). On each bootstrap sample the function 
    \code{\link{mledist}} is used to estimate bootstrapped values of parameters. When \code{\link{mledist}} fails
    to converge, \code{NA} values are returned. Medians and 2.5 and 97.5 percentiles are computed by removing
    \code{NA} values. The medians and the 95 percent confidence intervals of parameters (2.5 and 97.5 percentiles) 
    are printed in the summary.
    If inferior to the whole number of iterations, the number of iterations for which \code{\link{mledist}} converges 
    is also printed in the summary.
    
    The plot of an object of class \code{"bootdistcens"} consists in a scatterplot or a matrix of scatterplots
    of the bootstrapped values of parameters.
    It uses the function \code{\link{stripchart}} when the fitted distribution
    is characterized by only one parameter, and the function \code{\link{plot}} in other cases.
    In these last cases, it provides
    a representation of the joint uncertainty distribution of the fitted parameters.
    
        It is possible to accelerate the bootstrap using parallelization. We recommend you to
    use \code{parallel = "multicore"}, or \code{parallel = "snow"} if you work on Windows,
    and to fix \code{ncpus} to the number of available processors.

}

\value{ 
    \code{bootdistcens} returns an object of class \code{"bootdistcens"}, a list with 6 components,
    \item{estim}{ a data frame containing the bootstrapped values of parameters.}
    \item{converg}{ a vector containing the codes for convergence of the iterative method
        used to estimate parameters on 
        each bootstraped data set.}
    \item{method}{ A character string coding for the type of resampling : 
     in this case \code{"nonparam"} as it is the only available method for censored data. }
    \item{nbboot}{ The number of samples drawn by bootstrap.}   
    \item{CI}{ bootstrap medians and 95 percent confidence percentile intervals of parameters. }
    \item{fitpart}{ The object of class \code{"fitdistcens"} on which the bootstrap procedure
        was applied.}
    
Generic functions:
\describe{       
    \item{\code{print}}{
    The print of a \code{"bootdistcens"} object shows the bootstrap parameter estimates. If inferior to the whole number of bootstrap iterations, 
    the number of iterations 
    for which the estimation converges is also printed.
    }
    \item{\code{summary}}{
    The summary provides the median and 2.5 and 97.5 percentiles of each parameter. If inferior to the whole number of bootstrap iterations, 
    the number of iterations 
    for which the estimation converges is also printed in the summary.
    }
    \item{\code{plot}}{
    The plot shows the bootstrap estimates with the \code{\link{stripchart}} function
    for univariate parameters and \code{\link{plot}} function for multivariate parameters.
    }          
    
}   
}
        

\seealso{ 
    See \code{\link{fitdistrplus}} for an overview of the package.
    \code{\link{fitdistcens}}, \code{\link{mledist}}, \code{\link{quantile.bootdistcens}} 
    for another generic function to calculate 
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

%\note{  }

\examples{
# We choose a low number of bootstrap replicates in order to satisfy CRAN running times
# constraint.
# For practical applications, we recommend to use at least niter=501 or niter=1001.

# (1) Fit of a normal distribution to fluazinam data in log10
# followed by nonparametric bootstrap and calculation of quantiles
# with 95 percent confidence intervals
#
data(fluazinam)
(d1 <-log10(fluazinam))
f1 <- fitdistcens(d1, "norm")
b1 <- bootdistcens(f1, niter = 101)
b1
summary(b1)
plot(b1)
quantile(b1)
CIcdfplot(b1, CI.output = "quantile")

# (2) Estimation of the mean of the normal distribution 
# by maximum likelihood with the standard deviation fixed at 1 
# using the argument fix.arg
# followed by nonparametric bootstrap 
# and calculation of quantiles with 95 percent confidence intervals
#
f1b <- fitdistcens(d1, "norm", start = list(mean = 1),fix.arg = list(sd = 1))
b1b <- bootdistcens(f1b, niter = 101)
summary(b1b)
plot(b1b)
quantile(b1b)

# (3) comparison of sequential and parallel versions of bootstrap
# to be tried with a greater number of iterations (1001 or more)
#
\dontrun{
#small number to reduce computation time
niter <- 10 #in practice, choose above 100
data(fluazinam)
d1 <-log10(fluazinam)
f1 <- fitdistcens(d1, "norm")

# sequential version
ptm <- proc.time()
summary(bootdistcens(f1, niter = niter))
proc.time() - ptm

# parallel version using snow
require(parallel)
ptm <- proc.time()
summary(bootdistcens(f1, niter = niter, parallel = "snow", ncpus = 4))
proc.time() - ptm

# parallel version using multicore (not available on Windows)
ptm <- proc.time()
summary(bootdistcens(f1, niter = niter, parallel = "multicore", ncpus = 4))
proc.time() - ptm
}

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
