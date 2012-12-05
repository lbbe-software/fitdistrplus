\name{bootdist}
\alias{bootdist}
\alias{plot.bootdist}
\alias{print.bootdist}
\alias{summary.bootdist}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Bootstrap simulation of uncertainty for non-censored data}

\description{
 Uses parametric or nonparametric bootstrap resampling in order to simulate uncertainty 
 in the parameters of the distribution fitted to non-censored data.  
}

\usage{
bootdist(f, bootmethod="param", niter=1001)
\method{print}{bootdist}(x,...)
\method{plot}{bootdist}(x,...)
\method{summary}{bootdist}(object,...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{f}{ An object of class 'fitdist' result of the function \code{fitdist}.}
\item{bootmethod}{ A character string coding for the type of resampling : \code{"param"} for a parametric
    resampling and \code{"nonparam"} for a nonparametric resampling of data.}
\item{niter}{ The number of samples drawn by bootstrap.  }
\item{x}{ an object of class 'bootdist'.  }
\item{object}{ an object of class 'bootdist'.  }
\item{\dots}{ further arguments to be passed to generic methods }
}

\details{
    Samples are drawn by parametric bootstrap (resampling from the distribution fitted by \code{fitdist}) or
    non parametric bootstrap (resampling with replacement from the data set). On each bootstrap sample the function 
    \code{mledist}  (or \code{mmedist}, \code{qmedist} or \code{mgedist} 
    according to the component \code{f$method} of the 
    object of class 'fitdist') is used to estimate bootstrapped values of parameters. When that function fails
    to converge, \code{NA} values are returned. Medians and 2.5 and 97.5 percentiles are computed by removing
    \code{NA} values. 
    The medians and the 95 percent confidence intervals of parameters (2.5 and 97.5 percentiles) are printed 
    in the summary.
    If inferior to the whole number of iterations, the number of iterations for which the function converges 
    is also printed in the summary.
    
    The plot of an object of class \code{bootdist} consists in a scatterplot or a matrix of scatterplots
    of the bootstrapped values of parameters.
    It uses the function \code{\link{stripchart}} when the fitted distribution
    is characterized by only one parameter, and the function \code{\link{plot}} in other cases.
    In these last cases, it provides
    a representation of the joint uncertainty distribution of the fitted parameters.
}

\value{ 
    \code{bootdist} returns an object of class 'bootdist', a list with 4 components,
    \item{ estim }{ a data frame containing the bootstrapped values of parameters.}
    \item{ converg }{ a vector containing the codes for convergence obtained if an iterative method is
        used to estimate parameters on 
        each bootstraped data set (and 0 if a closed formula is used).}
    \item{ method }{ A character string coding for the type of resampling : \code{"param"} for a parametric
        resampling and \code{"nonparam"} for a nonparametric resampling of data. }
    \item{ CI }{ bootstrap medians and 95 percent confidence percentile intervals of parameters. }
}
        

\seealso{ 
    \code{\link{fitdist}}, \code{\link{mledist}}, \code{\link{qmedist}}, \code{\link{mmedist}}, and \code{\link{mgedist}}.
}

\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 181-241.
}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} 
}

%\note{  }

\examples{
# (1) Fit of a gamma distribution to serving size data
# using default method (maximum likelihood estimation)
# followed by parametric bootstrap
#
data(groundbeef)
x1 <- groundbeef$serving
f1 <- fitdist(x1, "gamma")
b1 <- bootdist(f1)
print(b1)
plot(b1)
summary(b1)

# (2) non parametric bootstrap on the same fit
# with less iterations
#
b1b <- bootdist(f1,bootmethod="nonparam",niter=501)
summary(b1b)

}


\keyword{ distribution }% at least one, from doc/KEYWORDS
