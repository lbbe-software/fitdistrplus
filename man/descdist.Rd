\name{descdist}
\alias{descdist}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Description of an empirical distribution for non-censored data}

\description{
 Computes descriptive parameters of an empirical distribution for non-censored data
 and provides a skewness-kurtosis plot.
}

\usage{
descdist(data,discrete=FALSE,boot=NULL,graph=TRUE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{data}{ A numeric vector.  }
\item{discrete}{ If \code{TRUE}, the distribution is considered as discrete.  }
\item{boot}{ If not \code{NULL}, \code{boot} values of skewness and kurtosis are plotted from bootstrap
 samples of \code{data}. \code{boot} must be fixed in this case to an integer above 10. }
\item{graph}{ If \code{FALSE}, the graph is not plotted. }
}

\details{
    Minimum, maximum, median, mean, sample sd, sample skewness and sample kurtosis values are printed. 
    A skewness-kurtosis plot such as the one proposed by Cullen and Frey (1999) is given for the 
    empirical distribution. On this plot, values for common distributions are also displayed as a tools 
    to help the choice of distributions to fit to data. In order to take into account the uncertainty 
    of the estimated values of kurtosis and skewness, the data set may be boostraped by 
    fixing the argument \code{boot} to an integer above 10. \code{boot} values of skewness and kurtosis 
    corresponding to the \code{boot} bootstrap samples are then computed and reported in blue color on the 
    skewness-kurtosis plot.
    
    If \code{discrete} is \code{TRUE}, 
    these common distributions are the Poisson, negative binomial and normal distributions. 
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
    \item{ sd }{ the standard deviation sample value  }
    \item{ skewness }{ the skewness sample value  }
    \item{ kurtosis }{ the kurtosis sample value  }    
}

\seealso{ 
    \code{\link{plotdist}} 
}

\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 81-159.
Evans M, Hastings N and Peacock B (2000) Statistical distributions. John Wiley and Sons Inc.
}

\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vet-lyon.fr} 
}

%\note{  }

\examples{
x1<-c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
descdist(x1)
descdist(x1,boot=1000)

x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
descdist(x2,discrete=TRUE)

x3<-rbeta(100,shape1=0.05,shape2=1)
descdist(x3,boot=1000)
}

\keyword{ distribution }% at least one, from doc/KEYWORDS
