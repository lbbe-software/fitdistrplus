\name{descdist}
\alias{descdist}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Description of an empirical distribution for non-censored data}

\description{
 Computes descriptive parameters of an empirical distribution for non-censored data
 and provides a skewness-kurtosis plot.
}

\usage{
descdist(data,discrete=FALSE,boot=NULL,method="unbiased",graph=TRUE,obs.col="red",boot.col="pink")
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{data}{ A numeric vector.  }
\item{discrete}{ If \code{TRUE}, the distribution is considered as discrete.  }
\item{boot}{ If not \code{NULL}, \code{boot} values of skewness and kurtosis are plotted from bootstrap
 samples of \code{data}. \code{boot} must be fixed in this case to an integer above 10. }
\item{method}{ "unbiased" for unbiased estimated values of statistics or "sample" for sample values. }
\item{graph}{ If \code{FALSE}, the skewness-kurtosis graph is not plotted. }
\item{obs.col}{ Color used for the observed point on the skewness-kurtosis graph. }
\item{boot.col}{ Color used for bootstrap sample of points on the skewness-kurtosis graph. }
}

\details{
    Minimum, maximum, median, mean, sample sd, and sample (if \code{method=="sample"}) or by default 
    unbiased estimations of skewness and 
    Pearsons's kurtosis values (Fisher, 1930) are printed. Be careful, estimations of skewness and kurtosis
    are unbiased only for normal distributions and estimated values are then only indicative.
    A skewness-kurtosis plot such as the one proposed by Cullen and Frey (1999) is given for the 
    empirical distribution. On this plot, values for common distributions are also displayed as a tools 
    to help the choice of distributions to fit to data. In order to take into account the uncertainty 
    of the estimated values of kurtosis and skewness, the data set may be boostraped by 
    fixing the argument \code{boot} to an integer above 10. \code{boot} values of skewness and kurtosis 
    corresponding to the \code{boot} bootstrap samples are then computed and reported in blue color on the 
    skewness-kurtosis plot.
    
    If \code{discrete} is \code{TRUE}, 
    the represented distributions are the Poisson, negative binomial and normal distributions. 
    If \code{discrete} is \code{FALSE}, these are uniform, normal, lognormal, beta and gamma distributions.
    The Weibull distribution is not represented on the graph but it is indicated on the legend that
    shapes close to lognormal and gamma distributions may be obtained with this distribution.
}

\value{ 
    \code{descdist} returns a list with 7 components,
    \item{ min }{ the minimum value  }
    \item{ max }{ the maximum value  }
    \item{ median }{ the median value  }
    \item{ mean }{ the mean value  }
    \item{ sd }{ the standard deviation sample or estimated value  }
    \item{ skewness }{ the skewness sample or estimated value  }
    \item{ kurtosis }{ the kurtosis sample or estimated value  }    
}

\seealso{ 
    \code{\link{plotdist}} 
}

\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 81-159.
Evans M, Hastings N and Peacock B (2000) Statistical distributions. John Wiley and Sons Inc.

Fisher RA (1930) The moments of the distribution for normal samples of measures of departures from normality.
Proc. R. Soc. London, Series A 130, 16-28.
}

\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vet-lyon.fr} 
}

%\note{  }

\examples{
x1 <- rnorm(100)
descdist(x1)
descdist(x1,boot=1000)

descdist(rbeta(100,shape1=0.05,shape2=1),boot=1000)
descdist(rgamma(100,shape=2,rate=1),boot=1000)

descdist(rpois(100,lambda=2),discrete=TRUE,boot=1000)
}

\keyword{ distribution }% at least one, from doc/KEYWORDS
