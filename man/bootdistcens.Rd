\name{bootdistcens}
\alias{bootdistcens}
\alias{plot.bootdistcens}
\alias{print.bootdistcens}
\alias{summary.bootdistcens}
\alias{quantile.bootdistcens}

\title{ Bootstrap simulation of uncertainty for censored data}

\description{
 Uses nonparametric bootstrap resampling in order to simulate uncertainty 
 in the parameters of the distribution fitted to censored data.  
}

\usage{
bootdistcens(f, niter=1001)
\method{print}{bootdistcens}(x,...)
\method{plot}{bootdistcens}(x,...)
\method{summary}{bootdistcens}(object,...)
\method{quantile}{bootdistcens}(x, probs = seq(0.1, 0.9, by=0.1), 
	CI.type = "two.sided", CI.level = 0.95, \dots)

}

\arguments{
\item{f}{ An object of class 'fitdistcens' result of the function \code{fitdistcens}.}
\item{niter}{ The number of samples drawn by bootstrap.  }
\item{x}{ An object of class 'bootdistcens'.  }
\item{object}{ An object of class 'bootdistcens'.  }
\item{\dots}{ Further arguments to be passed to generic methods }
\item{probs}{ A numeric vector of probabilities with values in [0, 1].}
\item{CI.type}{ Type of confidence interval: either \code{"two.sided"} or one-sided 
	intervals (\code{"less"} or \code{"greater"}).}
\item{CI.level}{ The confidence level required.}
}

\details{
    Samples are drawn by 
    non parametric bootstrap (resampling with replacement from the data set). On each bootstrap sample the function 
    \code{mledist} is used to estimate bootstrapped values of parameters. When \code{mledist} fails
    to converge, \code{NA} values are returned. Medians and 2.5 and 97.5 percentiles are computed by removing
    \code{NA} values. The medians and the 95 percent confidence intervals of parameters (2.5 and 97.5 percentiles) 
    are printed in the summary.
    If inferior to the whole number of iterations, the number of iterations for which \code{mledist} converges 
    is also printed in the summary.
    
    The plot of an object of class 'bootdistcens' consists in a scatterplot or a matrix of scatterplots
    of the bootstrapped values of parameters.
    It uses the function \code{\link{stripchart}} when the fitted distribution
    is characterized by only one parameter, and the function \code{\link{plot}} in other cases.
    In these last cases, it provides
    a representation of the joint uncertainty distribution of the fitted parameters.
}

\value{ 
    \code{bootdistcens} returns an object of class 'bootdistcens', a list with 3 components,
    \item{ estim }{ a data frame containing the bootstrapped values of parameters.}
    \item{ converg }{ a vector containing the codes for convergence obtained when using \code{mledist} on 
        each bootstraped data set.}
    \item{ CI }{ bootstrap medians and 95 percent confidence percentile intervals of parameters. }
	
Generic functions:
\describe{       
    \item{\code{print}}{
    The print of a 'bootdistcens' object shows the bootstrap estimates.
    }
    \item{\code{summary}}{
    The summary provides the median and the quartiles.
	}
    \item{\code{plot}}{
    The plot shows the bootstrap estimates with \code{stripchart} function
	for univariate parameters and \code{plot} function for multivariate parameters.
    }
    \item{\code{quantile}}{
    Quantiles of the parametric distribution are calculated for 
	each probability specified in \code{probs}, using the estimated parameters of \code{fitpart}.
    If \code{CI.type} is \code{two.sided}, 
    the \code{CI.level} two-sided confidence intervals of quantiles are calculated.
    If \code{CI.type} is \code{less}, 
    the \code{CI.level} one-sided confidence intervals of quantiles are calculated.
    If \code{CI.type} is \code{greater}, 
    the \code{CI.level} one-sided confidence intervals of quantiles are calculated.
    If inferior to the whole number of bootstrap iterations, the number of iterations 
	for which the estimation converges is also printed in the summary.
    }            
    
}	
}
        

\seealso{ 
    \code{\link{fitdistcens}} and \code{\link{mledist}}.
}

\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 181-241.
}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} 
}

%\note{  }

\examples{
# (1) Fit of a normal distribution to fluazinam data in log10
# followed by nonparametric bootstrap
#
data(fluazinam)
(d1 <-log10(fluazinam))
f1 <- fitdistcens(d1, "norm")
b1 <- bootdistcens(f1)
b1
summary(b1)
plot(b1)
quantile(b1)

# (2) Estimation of the mean of the normal distribution 
# by maximum likelihood with the standard deviation fixed at 1 
# using the argument fix.arg
# followed by nonparametric bootstrap with less iterations
#
f1b <- fitdistcens(d1, "norm", start = list(mean = 1),fix.arg = list(sd = 1))
b1b <- bootdistcens(f1b,niter = 501)
summary(b1b)
plot(b1b)
quantile(b1b)
quantile(b1b, probs = c(0.05, 0.1, 0.2), CI.type = "greater")
quantile(b1b, probs = c(0.05, 0.1, 0.2), CI.level = 0.90)

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
