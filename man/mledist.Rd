\name{mledist}
\alias{mledist}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Maximum likelihood fitting of univariate distributions}

\description{
 Fits a univariate distribution by maximum likelihood.  
}

\usage{
mledist(data, distr, start, optim.method="default", lower=-Inf, upper=Inf, optim.function=NULL, ...)
}

\arguments{
\item{data}{ A numeric vector. }
\item{distr}{ A character string \code{"name"} naming a distribution for which the corresponding 
    density function \code{dname} must be defined, or directly the density function. The first argument
	is assumed to be the argument for the vector of quantiles, such as in \code{pexp}.}
\item{start}{ A named list giving the initial values of parameters of the named distribution. 
    This argument may be omitted for some distributions for which reasonable 
    starting values are computed (see details).  }
\item{optim.method}{ \code{"default"} (see details) or optimization method to pass to \code{\link{optim}}.  }
\item{lower}{ Left bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{upper}{ Right bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{optim.function}{a function carrying the MLE optimisation, which must have a \code{fn} argument.}
\item{\dots}{further arguments passed to the \code{optim} or \code{optim.function} function.}
}

\details{
    When \code{optim.function=NULL} (the default), maximum likelihood estimations 
	of the distribution parameters are computed with the R base \code{\link{optim}}.
    Direct optimization of the log-likelihood is performed (using \code{\link{optim}}) by default with 
    the \code{"Nelder-Mead"} method for distributions characterized by more than one parameter and
    the \code{"BFGS"} method for distributions characterized by only one parameter, or 
    with the method specified in the argument \code{"optim.method"} if not \code{"default"}.
    Box-constrainted optimization may be used with the method \code{"L-BFGS-B"}, using the constraints
    on parameters specified in arguments \code{lower} and \code{upper}. If non-trivial bounds are supplied,
    this method will be automatically selected, with a warning.
    
    For the following named distributions, reasonable starting values will 
    be computed if \code{start} is omitted : \code{"norm"}, \code{"lnorm"},
    \code{"exp"} and \code{"pois"}, \code{"cauchy"}, \code{"gamma"}, \code{"logis"},
    \code{"nbinom"} (parametrized by mu and size), \code{"geom"}, \code{"beta"} and \code{"weibull"}. 
    Note that these starting 
    values may not be good enough if the fit is poor. The function is not able to 
	fit a uniform distribution.
    
    This function is not intended to be called directly but is internally called in
    \code{\link{fitdist}} and \code{\link{bootdist}} when used with the maximum likelihood method. 
}

\value{ 
    \code{mledist} returns a list with following components,
    \item{ estimate }{ the parameter estimates }
    \item{ convergence }{ an integer code for the convergence of \code{optim}. 
    
        \code{0} indicates successful convergence.
        
        \code{1} indicates that the iteration limit of \code{optim} has been reached.

        \code{10} indicates degeneracy of the Nealder-Mead simplex.

        \code{100} indicates that \code{optim} encountered an internal error.
        }
        
    \item{ loglik }{ the log-likelihood }
    \item{ hessian }{ a symmetric matrix computed by \code{optim} as an estimate of the Hessian
        at the solution found. It is used in \code{fitdist} to estimate standard errors. }
}


\seealso{ 
    \code{\link{momdist}}, \code{\link{fitdist}}, \code{\link{optim}} and \code{\link{bootdist}}.
}

\references{ 
    Venables WN and Ripley BD (2002) Modern applied statistics with S. Springer, New York, pp. 435-446.
}


\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vet-lyon.fr} 
}

\examples{
x1<-c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
mledist(x1,"norm")

dgumbel<-function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
mledist(x1,"gumbel",start=list(a=10,b=5))

x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
mledist(x2,"pois")
mledist(x2,"nbinom")

x3<-c(0.80,0.72,0.88,0.84,0.38,0.64,0.69,0.48,0.73,0.58,0.81,
0.83,0.71,0.75,0.59)
mledist(x3,"beta")

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
