\name{descdist}
\alias{descdist}
\alias{print.descdist}

\title{ Description of an empirical distribution for non-censored data}

\description{
 Computes descriptive parameters of an empirical distribution for non-censored data
 and provides a skewness-kurtosis plot.
}

\usage{
descdist(data, discrete = FALSE, boot = NULL, method = "unbiased",
graph = TRUE, print = TRUE, obs.col = "red", obs.pch = 16, boot.col = "orange")

\method{print}{descdist}(x, \dots)
}

\arguments{
\item{data}{ A numeric vector.  }
\item{discrete}{ If \code{TRUE}, the distribution is considered as discrete.  }
\item{boot}{ If not \code{NULL}, \code{boot} values of skewness and kurtosis are plotted from bootstrap
 samples of \code{data}. \code{boot} must be fixed in this case to an integer above 10. }
\item{method}{ "unbiased" for unbiased estimated values of statistics or "sample" for sample values. }
\item{graph}{ If \code{FALSE}, the skewness-kurtosis graph is not plotted. }
\item{print}{ If \code{FALSE}, the descriptive parameters computed are not printed. }
\item{obs.col}{ Color used for the observed point on the skewness-kurtosis graph. }
\item{obs.pch}{ plotting character used for the observed point on the skewness-kurtosis graph. }
\item{boot.col}{ Color used for bootstrap sample of points on the skewness-kurtosis graph. }
\item{x}{An object of class \code{"descdist"}.  }
\item{\dots}{Further arguments to be passed to generic functions}
}

\details{
    Minimum, maximum, median, mean, sample sd, and sample (if \code{method=="sample"}) or by default 
    unbiased estimations of skewness and 
    Pearsons's kurtosis values are printed (Sokal and Rohlf, 1995).
    A skewness-kurtosis plot such as the one proposed by Cullen and Frey (1999) is given for the 
    empirical distribution. On this plot, values for common distributions are also displayed as a tools 
    to help the choice of distributions to fit to data. For some distributions (normal, uniform,
    logistic, exponential for example), there is only one possible value for the skewness and the kurtosis
    (for a normal distribution for example, skewness = 0 and kurtosis = 3), and the distribution 
    is thus represented by a point on the plot. For other distributions, 
    areas of possible values are represented, consisting in lines (gamma and lognormal distributions for example), 
    or larger areas (beta distribution for example). The Weibull distribution is not represented on the graph but it 
    is indicated on the legend that
    shapes close to lognormal and gamma distributions may be obtained with this distribution.
    
    In order to take into account the uncertainty 
    of the estimated values of kurtosis and skewness from data, the data set may be bootstraped by 
    fixing the argument \code{boot} to an integer above 10. \code{boot} values of skewness and kurtosis 
    corresponding to the \code{boot} bootstrap samples are then computed and reported in blue color on the 
    skewness-kurtosis plot.
    
    If \code{discrete} is \code{TRUE}, 
    the represented distributions are the Poisson, negative binomial distributions, 
    and the normal distribution to which previous discrete distributions may converge. 
    If \code{discrete} is \code{FALSE}, these are uniform, normal, logistic, lognormal, beta 
    and gamma distributions.
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
    \item{method}{ the method specified in input ("unbiased" for unbiased 
    estimated values of statistics or "sample" for sample values. }
}

\seealso{ 
    See \code{\link{plotdist}}.
    
    Please visit the \href{https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html}{Frequently Asked Questions}.
}

\references{ 
Cullen AC and Frey HC (1999), \emph{Probabilistic techniques in exposure assessment}.
Plenum Press, USA, pp. 81-159.

Evans M, Hastings N and Peacock B (2000), \emph{Statistical distributions}.
John Wiley and Sons Inc, \doi{10.1002/9780470627242}.

Sokal RR and Rohlf FJ (1995), \emph{Biometry}.
W.H. Freeman and Company, USA, pp. 111-115.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34, \doi{https://doi.org/10.18637/jss.v064.i04}.

}

\author{ 
Marie-Laure Delignette-Muller and Christophe Dutang.
}

\examples{
# (1) Description of a sample from a normal distribution
# with and without uncertainty on skewness and kurtosis estimated by bootstrap 
#
set.seed(1234)
x1 <- rnorm(100)
descdist(x1)
descdist(x1,boot=11)

# (2) Description of a sample from a beta distribution
# with uncertainty on skewness and kurtosis estimated by bootstrap
# with changing of default colors and plotting character for observed point
#
descdist(rbeta(100,shape1=0.05,shape2=1),boot=11,
 obs.col="blue", obs.pch = 15, boot.col="darkgreen")

# (3) Description of a sample from a gamma distribution
# with uncertainty on skewness and kurtosis estimated by bootstrap
# without plotting 
#
descdist(rgamma(100,shape=2,rate=1),boot=11,graph=FALSE)

# (4) Description of a sample from a Poisson distribution
# with uncertainty on skewness and kurtosis estimated by bootstrap 
#
descdist(rpois(100,lambda=2),discrete=TRUE,boot=11)

# (5) Description of serving size data
# with uncertainty on skewness and kurtosis estimated by bootstrap 
#
data(groundbeef)
serving <- groundbeef$serving
descdist(serving, boot=11)
}

\keyword{ distribution }
