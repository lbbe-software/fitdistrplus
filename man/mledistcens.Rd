\name{mledistcens}
\alias{mledistcens}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Maximum likelihood fitting of univariate distributions from censored data}

\description{
 Fits a univariate distribution from censored data 
 by maximum likelihood using the optimization function \code{\link{optim}}.  
}

\usage{
mledistcens(censdata, distr, start)
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
}

\details{
    Maximum likelihood estimations of the distribution parameters are computed. 
    Direct optimization of the log-likelihood is performed using \code{\link{optim}}, with its
    default method \code{"Nelder-Mead"} for distributions characterized by more than one parameter and
    the method \code{"BFGS"} for distributions characterized by only one parameter. 
    For the following named distributions, reasonable starting values will 
    be computed if \code{start} is omitted : \code{"norm"}, \code{"lnorm"},
    \code{"exp"} and \code{"pois"}, \code{"cauchy"}, \code{"gamma"}, \code{"logis"},
    \code{"nbinom"} (parametrized by mu and size), \code{"geom"}, \code{"beta"} and \code{"weibull"}. 
    Note that these starting 
    values may not be good enough if the fit is poor. The function is not able to fit a uniform distribution.
    
    This function is not intended to be called directly but is internally called in
    \code{\link{fitdistcens}} and \code{\link{bootdistcens}}. 
}

\value{ 
    \code{mledistcens} returns a list with 4 components,
    \item{ estimate }{ the parameter estimates }
    \item{ convergence }{ an integer code for the convergence of \code{optim}. 
    
        \code{0} indicates successful convergence.
        
        \code{1} indicates that the iteration limit of \code{optim} has been reached.

        \code{10} indicates degeneracy of the Nealder-Mead simplex.

        \code{100} indicates that \code{optim} encountered an internal error.
        }
        
    \item{ loglik }{ the log-likelihood }
    \item{ hessian }{ a symmetric matrix computed by \code{optim} as an estimate of the Hessian
        at the solution found. It is used in \code{fitdistcens} to estimate standard errors. }
}


\seealso{ 
    \code{\link{fitdistcens}}, \code{\link{optim}} and \code{\link{bootdistcens}}.
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
left=c(1.73,1.51,0.77,1.96,1.96,-1.4,-1.4,NA,-0.11,0.55,0.41,2.56,NA,-0.53,0.63,-1.4,-1.4,-1.4,NA,0.13),
right=c(1.73,1.51,0.77,1.96,1.96,0,-0.7,-1.4,-0.11,0.55,0.41,2.56,-1.4,-0.53,0.63,0,-0.7,NA,-1.4,0.13))
mledistcens(d1,"norm")

dgumbel<-function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel<-function(q,a,b) exp(-exp((a-q)/b))
mledistcens(d1,"gumbel",start=list(a=0,b=2))

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
