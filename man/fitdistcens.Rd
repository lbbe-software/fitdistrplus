\name{fitdistcens}
\alias{fitdistcens}
\alias{plot.fitdistcens}
\alias{print.fitdistcens}
\alias{summary.fitdistcens}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fitting of univariate distributions to censored data}

\description{
 Fits a univariate distribution to censored data by maximum likelihood.  
}

\usage{
fitdistcens(censdata, distr, start,...)
\method{print}{fitdistcens}(x,...)
\method{plot}{fitdistcens}(x,...)
\method{summary}{fitdistcens}(object,...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{censdata}{ A dataframe of two columns respectively named \code{left} 
    and \code{right}, describing each observed value as an interval.
    The \code{left} column contains either \code{NA} for left censored observations,
     the left bound of the interval for interval censored observations,
     or the observed value for non-censored observations.
     The \code{right} column contains either \code{NA} for right censored observations,
     the right bound of the interval for interval censored observations,
     or the observed value for non-censored observations.  }
\item{distr}{ A character string \code{"name"} naming a distribution, for which the corresponding 
    density function \code{dname} and the corresponding distribution function \code{pname}
    must be defined, or directly the density function.}
\item{start}{ A named list giving the initial values of parameters of the named distribution. 
    This argument may be omitted for some distributions for which reasonable 
    starting values are computed (see details).  }
\item{x}{ an object of class 'fitdistcens'.  }
\item{object}{ an object of class 'fitdistcens'.  }
\item{...}{ further arguments to be passed to generic functions, or to the function \code{"mledist"} 
     in order to control the optimization method.}
}

\details{
    Maximum likelihood estimations of the distribution parameters are computed using 
    the function \code{\link{mledist}}. 
    By default direct optimization of the log-likelihood is performed using \code{\link{optim}},
    with the "Nelder-Mead" method for distributions characterized by more than one parameter
    and the "BFGS" method for distributions characterized by only one parameter. 
    The method used in \code{\link{optim}} may be chosen or another optimization method
    may be chosen using ... argument (see \code{\link{mledist}} for details).
    For the following named distributions, reasonable starting values will 
    be computed if \code{start} is omitted : \code{"norm"}, \code{"lnorm"},
    \code{"exp"} and \code{"pois"}, \code{"cauchy"}, \code{"gamma"}, \code{"logis"},
    \code{"nbinom"} (parametrized by mu and size), \code{"geom"}, \code{"beta"} and \code{"weibull"}. 
    Note that these starting 
    values may not be good enough if the fit is poor. The function is not able to fit a uniform distribution.     
    With the parameter estimates, the function returns the log-likelihood and the standard errors of 
    the estimates calculated from the 
    Hessian at the solution found by \code{optim}.
     
    The plot of an object of class "fitdistcens" returned by \code{fitdistcens} uses the function \code{\link{plotdistcens}}.
}

\value{ 
    \code{fitdistcens} returns an object of class 'fitdistcens', a list with following components,
    \item{ estimate }{ the parameter estimates }
    \item{ sd }{ the estimated standard errors }
    \item{ cor }{ the estimated correlation matrix }
    \item{ loglik }{ the log-likelihood }
     \item{ aic }{ the Akaike information criterion }
    \item{ bic }{ the the so-called BIC or SBC (Schwarz Bayesian criterion)}
    \item{ censdata }{ the censored dataset }
    \item{ distname }{ the name of the distribution }
}


\seealso{ 
    \code{\link{plotdistcens}}, \code{\link{optim}}, \code{\link{mledist}} and 
    \code{\link{fitdist}}.
}

\references{ 
Venables WN and Ripley BD (2002) Modern applied statistics with S. Springer, New York, pp. 435-446.
}


\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vet-lyon.fr} 
}

%\note{  }

\examples{
d1<-data.frame(
left=c(1.73,1.51,0.77,1.96,1.96,-1.4,-1.4,NA,-0.11,0.55,0.41,
    2.56,NA,-0.53,0.63,-1.4,-1.4,-1.4,NA,0.13),
right=c(1.73,1.51,0.77,1.96,1.96,0,-0.7,-1.4,-0.11,0.55,0.41,
    2.56,-1.4,-0.53,0.63,0,-0.7,NA,-1.4,0.13))
f1n<-fitdistcens(d1, "norm")
f1n
summary(f1n)
plot(f1n,rightNA=3)

dgumbel<-function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel<-function(q,a,b) exp(-exp((a-q)/b))
f1g<-fitdistcens(d1,"gumbel",start=list(a=0,b=2))
summary(f1g)
plot(f1g,rightNA=3)

d3<-data.frame(left=10^(d1$left),right=10^(d1$right))
f3w<-fitdistcens(d3,"weibull")
summary(f3w)
plot(f3w,leftNA=0)
f3l<-fitdistcens(d3,"lnorm")
summary(f3l)
plot(f3l,leftNA=0)
f3e<-fitdistcens(d3,"exp")
summary(f3e)
plot(f3e,leftNA=0)

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
