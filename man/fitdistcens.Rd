\name{fitdistcens}
\alias{fitdistcens}
\alias{plot.fitdistcens}
\alias{print.fitdistcens}
\alias{summary.fitdistcens}
\title{Fitting of univariate distributions to censored data}

\description{
 Fits a univariate distribution to censored data by maximum likelihood.  
}

\usage{
fitdistcens(censdata, distr, start=NULL, fix.arg=NULL, \dots)

\method{print}{fitdistcens}(x, \dots)

\method{plot}{fitdistcens}(x, \dots)

\method{summary}{fitdistcens}(object, \dots)

    
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
    that must be kept fixed rather than estimated by maximum likelihood. }
\item{x}{ an object of class \code{"fitdistcens"}.  }
\item{object}{ an object of class \code{"fitdistcens"}.  }
\item{\dots}{ further arguments to be passed to generic functions,
    to the function \code{plotdistcens} in order to control the type of ecdf-plot used for censored data, 
    or to the function \code{mledist} 
    in order to control the optimization method.}
}

\details{
    Maximum likelihood estimations of the distribution parameters are computed using 
    the function \code{\link{mledist}}. 
    By default direct optimization of the log-likelihood is performed using \code{\link{optim}},
    with the "Nelder-Mead" method for distributions characterized by more than one parameter
    and the "BFGS" method for distributions characterized by only one parameter. 
    The algorithm used in \code{\link{optim}} can be chosen or another optimization function
    can be specified using \dots argument (see \code{\link{mledist}} for details).
    For the following named distributions, reasonable starting values will 
    be computed if \code{start} is omitted : \code{"norm"}, \code{"lnorm"},
    \code{"exp"} and \code{"pois"}, \code{"cauchy"}, \code{"gamma"}, \code{"logis"},
    \code{"nbinom"} (parametrized by mu and size), \code{"geom"}, \code{"beta"} and \code{"weibull"}. 
    Note that these starting 
    values may not be good enough if the fit is poor. The function is not able to fit a uniform distribution.     
    With the parameter estimates, the function returns the log-likelihood and the standard errors of 
    the estimates calculated from the 
    Hessian at the solution found by \code{optim} or by the user-supplied function passed to mledist.
     

}

\value{ 
    \code{fitdistcens} returns an object of class \code{"fitdistcens"}, a list with the following components:
    \item{estimate}{ the parameter estimates }
    \item{sd}{ the estimated standard errors }
    \item{cor}{ the estimated correlation matrix }
    \item{loglik}{ the log-likelihood }
    \item{aic}{ the Akaike information criterion }
    \item{bic}{ the the so-called BIC or SBC (Schwarz Bayesian criterion)}
    \item{censdata}{ the censored data set }
    \item{distname}{ the name of the distribution }
    \item{fix.arg}{ the named list giving the values of parameters of the named distribution
    that must be kept fixed rather than estimated by maximum likelihood or NULL if there are no such parameters. }
    \item{dots}{ the list of  further arguments passed in \dots to be used in \code{bootdistcens} to
    control the optimization method used in iterative calls to \code{mledist} or NULL if no such arguments.}
    
    
    Generic functions:
\describe{       
    \item{\code{print}}{
    The print of a \code{"fitdist"} object shows few traces about the fitting method and the fitted distribution.
    }
    \item{\code{summary}}{
    The summary provides the parameter estimates of the fitted distribution, the log-likelihood, 
    AIC and BIC statistics, the standard errors of the 
    parameter estimates and the correlation matrix between parameter estimates.
    }
    \item{\code{plot}}{
    The plot of an object of class \code{"fitdistcens"} returned by \code{fitdistcens} uses the 
    function \code{\link{plotdistcens}}.
    }
}
}


\seealso{ 
    \code{\link{plotdistcens}}, \code{\link{optim}}, \code{\link{mledist}}, 
    \code{\link{fitdist}} and \code{\link{quantile.fitdistcens}} for another generic function to calculate 
    quantiles from the fitted distribution.
}

\references{ 
Venables WN and Ripley BD (2002) Modern applied statistics with S. Springer, New York, pp. 435-446.
}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr}.
}

%\note{  }

\examples{
# (1) Fit of a lognormal distribution to bacterial contamination data
#
data(smokedfish)
fitsf  <-  fitdistcens(smokedfish,"lnorm")
summary(fitsf)
# default plot using the Turnbull algorithm
plot(fitsf)
# plot using the Turnbull algorithm with confidence intervals for the 
# empirical distribution
plot(fitsf, Turnbull.confint = TRUE)
# basic plot using intervals and points (see ?plotdiscens for details)
plot(fitsf, Turnbull = FALSE)
# plot of the same fit using the Turnbull algorithm in logscale
cdfcompcens(fitsf,main="bacterial contamination fits",
    xlab="bacterial concentration (CFU/g)",ylab="F",
    addlegend = FALSE,lines01 = TRUE, xlogscale = TRUE, xlim = c(1e-2,1e2))
# zoom on large values of F
cdfcompcens(fitsf,main="bacterial contamination fits",
    xlab="bacterial concentration (CFU/g)",ylab="F",
    addlegend = FALSE,lines01 = TRUE, xlogscale = TRUE, 
    xlim = c(1e-2,1e2),ylim=c(0.4,1))

# (2) Fit of a normal distribution on acute toxicity values 
# of fluazinam (in decimal logarithm) for
# macroinvertebrates and zooplancton, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology
#

data(fluazinam)
log10EC50 <-log10(fluazinam)
fln <- fitdistcens(log10EC50,"norm")
fln
summary(fln)
plot(fln)

# (3) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view dedicated to 
# probability distributions
#

dgumbel  <-  function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel  <-  function(q,a,b) exp(-exp((a-q)/b))
qgumbel  <-  function(p,a,b) a-b*log(-log(p))
fg <- fitdistcens(log10EC50,"gumbel",start=list(a=1,b=1))
summary(fg)
plot(fg)

# (4) comparison of fits of various distributions
# 

fll <- fitdistcens(log10EC50,"logis")
summary(fll)

cdfcompcens(list(fln,fll,fg),legendtext=c("normal","logistic","gumbel"),
xlab = "log10(EC50)")

# (5) how to change the optimisation method?
#

fitdistcens(log10EC50,"logis",optim.method="Nelder-Mead")
fitdistcens(log10EC50,"logis",optim.method="BFGS") 
fitdistcens(log10EC50,"logis",optim.method="SANN") 

# (6) custom optimisation function - example with the genetic algorithm
#
\dontrun{

    #wrap genoud function rgenoud package
    mygenoud  <-  function(fn, par, ...) 
    {
        require(rgenoud)
        res  <-  genoud(fn, starting.values=par, ...)        
        standardres  <-  c(res, convergence=0)
            
        return(standardres)
    }

    # call fitdistcens with a 'custom' optimization function
    fit.with.genoud <- fitdistcens(log10EC50,"logis", custom.optim=mygenoud, nvars=2,    
        Domains=cbind(c(0,0), c(5, 5)), boundary.enforcement=1, 
        print.level=1, hessian=TRUE)

    summary(fit.with.genoud)
}

# (7) estimation of the mean of a normal distribution 
# by maximum likelihood with the standard deviation fixed at 1 using the argument fix.arg
#
flnb <- fitdistcens(log10EC50, "norm", start = list(mean = 1),fix.arg = list(sd = 1))

# (8) Fit of a lognormal distribution on acute toxicity values of fluazinam for
# macroinvertebrates and zooplancton, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5 percent quantile value of 
# the fitted distribution (which is called the 5 percent hazardous concentration, HC5,
# in ecotoxicology) and estimation of other quantiles.

data(fluazinam)
log10EC50 <-log10(fluazinam)
fln <- fitdistcens(log10EC50,"norm")

quantile(fln, probs = 0.05)
quantile(fln, probs = c(0.05, 0.1, 0.2))

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
