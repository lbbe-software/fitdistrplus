\name{fitdistcens}
\alias{fitdistcens}
\alias{plot.fitdistcens}
\alias{print.fitdistcens}
\alias{summary.fitdistcens}
\alias{logLik.fitdistcens}
\alias{vcov.fitdistcens}
\alias{coef.fitdistcens}
\title{Fitting of univariate distributions to censored data}

\description{
 Fits a univariate distribution to censored data by maximum likelihood.  
}

\usage{
fitdistcens(censdata, distr, start=NULL, fix.arg=NULL, 
  keepdata = TRUE, keepdata.nb=100, \dots)

\method{print}{fitdistcens}(x, \dots)

\method{plot}{fitdistcens}(x, \dots)

\method{summary}{fitdistcens}(object, \dots)

\method{logLik}{fitdistcens}(object, \dots)

\method{vcov}{fitdistcens}(object, \dots)

\method{coef}{fitdistcens}(object, \dots)

    
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
    starting values are computed (see the 'details' section of  \code{\link{mledist}}).  }
\item{fix.arg}{ An optional named list giving the values of parameters of the named distribution
    that must be kept fixed rather than estimated by maximum likelihood. }
\item{x}{ an object of class \code{"fitdistcens"}.  }
\item{object}{ an object of class \code{"fitdistcens"}.  }

\item{keepdata}{a logical. If \code{TRUE}, dataset is returned, 
otherwise only a sample subset is returned.}
\item{keepdata.nb}{When \code{keepdata=FALSE}, the length of the subset returned.}

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
    \code{start} may be omitted (i.e. \code{NULL}) for some classic distributions
    (see the 'details' section of  \code{\link{mledist}}).
        Note that when errors are raised by \code{optim}, it's a good idea to start by adding traces during
    the optimization process by adding \code{control=list(trace=1, REPORT=1)} in \dots argument.

    The function is not able to fit a uniform distribution.     
    With the parameter estimates, the function returns the log-likelihood and the standard errors of 
    the estimates calculated from the 
    Hessian at the solution found by \code{optim} or by the user-supplied function passed to mledist.
    
    By default (\code{keepdata = TRUE}), the object returned by \code{fitdist} contains 
    the data vector given in input.
    When dealing with large datasets, we can remove the original dataset from the output by
    setting \code{keepdata = FALSE}. In such a case, only \code{keepdata.nb} points (at most) 
    are kept by random subsampling \code{keepdata.nb}-4 points from the dataset and 
    adding the component-wise minimum and maximum. 
    If combined with \code{\link{bootdistcens}}, be aware that bootstrap is performed on the subset
    randomly selected in \code{fitdistcens}. Currently, the graphical comparisons of multiple fits 
    is not available in this framework.
    
    Weighted version of the estimation process is available for \code{method = "mle"}
    by using \code{weights=\dots}. See the corresponding man page for details.
    It is not yet possible to take into account weighths in functions plotdistcens, 
    plot.fitdistcens and cdfcompcens 
    (developments planned in the future).

}

\value{ 
    \code{fitdistcens} returns an object of class \code{"fitdistcens"}, a list with the following components:
    \item{estimate}{ the parameter estimates.}
    \item{method}{ the character string coding for the fitting method : 
        only \code{"mle"} for 'maximum likelihood estimation'.}
    \item{sd}{ the estimated standard errors.}
    \item{cor}{ the estimated correlation matrix, \code{NA} if numerically not computable
    or \code{NULL} if not available.}
    \item{vcov}{ the estimated variance-covariance matrix, \code{NULL} if not available.}
    \item{loglik}{ the log-likelihood.}
    \item{aic}{ the Akaike information criterion.}
    \item{bic}{ the the so-called BIC or SBC (Schwarz Bayesian criterion).}
    \item{censdata}{ the censored data set.}
    \item{distname}{ the name of the distribution.}
    \item{fix.arg}{ the named list giving the values of parameters of the named distribution
      that must be kept fixed rather than estimated by maximum likelihood or 
      \code{NULL} if there are no such parameters. }
    \item{fix.arg.fun}{the function used to set the value of \code{fix.arg} or \code{NULL}.}    
    \item{dots}{ the list of  further arguments passed in \dots to be used in \code{bootdistcens} to
    control the optimization method used in iterative calls to \code{mledist} or 
      \code{NULL} if no such arguments.}
    \item{convergence}{ an integer code for the convergence of 
    \code{\link{optim}}/\code{\link{constrOptim}} defined as below
    or defined by the user in the user-supplied optimization function. 
        \code{0} indicates successful convergence.
        \code{1} indicates that the iteration limit of \code{\link{optim}} has been reached.
        \code{10} indicates degeneracy of the Nealder-Mead simplex.
        \code{100} indicates that \code{\link{optim}} encountered an internal error.} 
    \item{discrete}{always \code{FALSE}.}    
    \item{weights}{the vector of weigths used in the estimation process or \code{NULL}.}
    

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
    \item{\code{logLik}}{ Extracts the estimated log-likelihood from the \code{"fitdistcens"} object.
    }
    \item{\code{vcov}}{ Extracts the estimated var-covariance matrix from the \code{"fitdistcens"} object
    (only available When \code{method = "mle"}).
    }
    \item{\code{coef}}{ Extracts the fitted coefficients from the \code{"fitdistcens"} object.
    }
}
}


\seealso{ 
    See \code{\link{fitdistrplus}} for an overview of the package.
    \code{\link{plotdistcens}}, \code{\link{optim}}, \code{\link{mledist}}, 
    \code{\link{fitdist}} and \code{\link{quantile.fitdistcens}} for another generic function to calculate 
    quantiles from the fitted distribution.
}

\references{ 
Venables WN and Ripley BD (2002), \emph{Modern applied statistics with S}.
Springer, New York, pp. 435-446.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}


\author{ 
Marie-Laure Delignette-Muller and Christophe Dutang.
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
