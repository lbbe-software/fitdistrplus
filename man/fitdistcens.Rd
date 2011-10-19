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
fitdistcens(censdata, distr, start=NULL, fix.arg=NULL,...)
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
\item{fix.arg}{ An optional named list giving the values of parameters of the named distribution
    that must kept fixed rather than estimated by maximum likelihood. }
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
    Hessian at the solution found by \code{optim} or by the user-supplied function passed to mledist.
     
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
    \item{ fix.arg }{ the named list giving the values of parameters of the named distribution
    that must kept fixed rather than estimated by maximum likelihood or NULL if there are no such parameters. }
    \item{ dots }{ the list of  further arguments passed in \dots to be used in \code{bootdistcens} to
    control the optimization method used in iterative calls to \code{mledist} or NULL if no such arguments}
}


\seealso{ 
    \code{\link{plotdistcens}}, \code{\link{optim}}, \code{\link{mledist}} and 
    \code{\link{fitdist}}.
}

\references{ 
Venables WN and Ripley BD (2002) Modern applied statistics with S. Springer, New York, pp. 435-446.
}


\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vetagro-sup.fr} 
}

%\note{  }

\examples{
# (1) basic fit of a normal distribution on censored data
#

d1<-data.frame(
left=c(1.73,1.51,0.77,1.96,1.96,-1.4,-1.4,NA,-0.11,0.55,0.41,
    2.56,NA,-0.53,0.63,-1.4,-1.4,-1.4,NA,0.13),
right=c(1.73,1.51,0.77,1.96,1.96,0,-0.7,-1.4,-0.11,0.55,0.41,
    2.56,-1.4,-0.53,0.63,0,-0.7,NA,-1.4,0.13))
f1n<-fitdistcens(d1, "norm")
f1n
summary(f1n)
plot(f1n)

# (2) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view dedicated to 
# probability distributions
#

dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q,a,b) exp(-exp((a-q)/b))
qgumbel <- function(p,a,b) a-b*log(-log(p))
f1g<-fitdistcens(d1,"gumbel",start=list(a=0,b=2))
summary(f1g)
plot(f1g)

# (3) comparison of fits of various distributions
# 

d3<-data.frame(left=10^(d1$left),right=10^(d1$right))
f3w<-fitdistcens(d3,"weibull")
summary(f3w)
plot(f3w)
f3l<-fitdistcens(d3,"lnorm")
summary(f3l)
plot(f3l)
f3e<-fitdistcens(d3,"exp")
summary(f3e)
plot(f3e)

# (4) how to change the optimisation method?
#

fitdistcens(d3,"gamma",optim.method="Nelder-Mead")
fitdistcens(d3,"gamma",optim.method="BFGS") 
fitdistcens(d3,"gamma",optim.method="SANN") 
fitdistcens(d3,"gamma",optim.method="L-BFGS-B",lower=c(0,0)) 

# (5) custom optimisation function - example with the genetic algorithm
#
\dontrun{

    #wrap genoud function rgenoud package
    mygenoud <- function(fn, par, ...) 
    {
        require(rgenoud)
        res <- genoud(fn, starting.values=par, ...)        
        standardres <- c(res, convergence=0)
            
        return(standardres)
    }

    # call fitdistcens with a 'custom' optimization function
    fit.with.genoud<-fitdistcens(d3, "gamma", custom.optim=mygenoud, nvars=2,    
        Domains=cbind(c(0,0), c(10, 10)), boundary.enforcement=1, 
        print.level=1, hessian=TRUE)

    summary(fit.with.genoud)
}

# (6) estimation of the standard deviation of a normal distribution 
# by maximum likelihood with the mean fixed at 0.1 using the argument fix.arg
#
fitdistcens(d1, "norm", start=list(sd=1.5),fix.arg=list(mean=0.1))

# (7) Fit of a lognormal distribution to bacterial contamination data
#
data(smokedfish)
fitsf <- fitdistcens(smokedfish,"norm")
summary(fitsf)
plot(fitsf)
plot(fitsf, Turnbull = FALSE)

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
