\name{mledist}
\alias{mledist}
\alias{mle}
\title{ Maximum likelihood fit of univariate distributions}

\description{
 Fit of univariate distributions using maximum likelihood for censored or non censored data.  
}

\usage{
mledist(data, distr, start = NULL, fix.arg = NULL, optim.method = "default", 
    lower = -Inf, upper = Inf, custom.optim = NULL, \dots)
}

\arguments{
\item{data}{ A numeric vector for non censored data or
    a dataframe of two columns respectively named \code{left} 
    and \code{right}, describing each observed value as an interval for censored data.
    In that case the \code{left} column contains either \code{NA} for left censored observations,
     the left bound of the interval for interval censored observations,
     or the observed value for non-censored observations.
     The \code{right} column contains either \code{NA} for right censored observations,
     the right bound of the interval for interval censored observations,
     or the observed value for non-censored observations. }
\item{distr}{ A character string \code{"name"} naming a distribution 
     for which the corresponding 
    density function \code{dname} and the corresponding distribution \code{pname}
    must be classically defined. }
\item{start}{ A named list giving the initial values of parameters of the named distribution. 
    This argument may be omitted for some distributions for which reasonable 
    starting values are computed (see details).  }
\item{fix.arg}{ An optional named list giving the values of parameters of the named distribution
    that must kept fixed rather than estimated by maximum likelihood. }
\item{optim.method}{ \code{"default"} (see details) or optimization method to pass to \code{\link{optim}}.  }
\item{lower}{ Left bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{upper}{ Right bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{custom.optim}{a function carrying the MLE optimisation (see details).}
\item{\dots}{further arguments passed to the \code{optim} or \code{custom.optim} function.}
}

\details{
    When \code{custom.optim=NULL} (the default), maximum likelihood estimations 
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
    
    If \code{custom.optim} is not \code{NULL}, then the user-supplied function is used instead
    of the R base \code{optim}. The \code{custom.optim} must have (at least) the following arguments
    \code{fn} for the function to be optimized, \code{par} for the initialized parameters. Internally the
    function to be optimized will also have other arguments, such as \code{obs} with observations and \code{ddistname}
    with distribution name for non censored data (Beware of potential conflicts with optional arguments 
    of \code{custom.optim}). It is
    assumed that \code{custom.optim} should carry out a MINIMIZATION. Finally, it should return
    at least the following components \code{par} for the estimate, \code{convergence} for the convergence
    code, \code{value} for \code{fn(par)} and \code{hessian}. See examples in \code{\link{fitdist}}
    and \code{\link{fitdistcens}}.
    
    This function is not intended to be called directly but is internally called in
    \code{\link{fitdist}} and \code{\link{bootdist}} when used with the maximum likelihood method
    and \code{\link{fitdistcens}} and \code{\link{bootdistcens}}. 
    
    NB: if your data values are particularly small or large, a scaling may be needed 
    before the optimization process. See example (7).
}

\value{ 
    \code{mledist} returns a list with following components,
    \item{estimate}{ the parameter estimates }
    \item{convergence}{ an integer code for the convergence of \code{optim} defined as below
    or defined by the user in the user-supplied optimization function. 
    
        \code{0} indicates successful convergence.
        
        \code{1} indicates that the iteration limit of \code{optim} has been reached.

        \code{10} indicates degeneracy of the Nealder-Mead simplex.

        \code{100} indicates that \code{optim} encountered an internal error.
        }
        
    \item{loglik}{ the log-likelihood }
    \item{hessian}{ a symmetric matrix computed by \code{optim} as an estimate of the Hessian
        at the solution found or computed in the user-supplied optimization function.
        It is used in \code{fitdist} to estimate standard errors. }
    \item{optim.function}{ the name of the optimization function used for maximum likelihood}
    \item{fix.arg}{ the named list giving the values of parameters of the named distribution
        that must kept fixed rather than estimated by maximum likelihood or \code{NULL} 
        if there are no such parameters. }

}


\seealso{ 
    \code{\link{mmedist}}, \code{\link{qmedist}}, \code{\link{mgedist}}, 
    \code{\link{fitdist}},\code{\link{fitdistcens}}, \code{\link{optim}},
    \code{\link{bootdistcens}} and \code{\link{bootdist}}.
}

\references{ 
    Venables W.N. and Ripley B.D. (2002) Modern applied statistics with S. Springer, New York, pp. 435-446.
}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} and
Christophe Dutang.
}

\examples{

# (1) basic fit of a normal distribution with maximum likelihood estimation
#

set.seed(1234)
x1 <- rnorm(n=100)
mledist(x1,"norm")

# (2) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view dedicated to probability distributions

dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
mledist(x1,"gumbel",start=list(a=10,b=5))

# (3) fit of a discrete distribution (Poisson)
#

set.seed(1234)
x2 <- rpois(n=30,lambda = 2)
mledist(x2,"pois")

# (4) fit a finite-support distribution (beta)
#

set.seed(1234)
x3 <- rbeta(n=100,shape1=5, shape2=10)
mledist(x3,"beta")


# (5) fit frequency distributions on USArrests dataset.
#

x4 <- USArrests$Assault
mledist(x4, "pois")
mledist(x4, "nbinom")

# (6) fit a continuous distribution (Gumbel) to censored data.
#

data(fluazinam)
log10EC50 <-log10(fluazinam)
# definition of the Gumbel distribution
dgumbel  <-  function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel  <-  function(q,a,b) exp(-exp((a-q)/b))
qgumbel  <-  function(p,a,b) a-b*log(-log(p))

mledist(log10EC50,"gumbel",start=list(a=0,b=2),optim.method="Nelder-Mead")

# (7) scaling problem
# the simulated dataset (below) has particularly small values, hence without scaling (10^0),
# the optimization raises an error. The for loop shows how scaling by 10^i
# for i=1,...,6 makes the fitting procedure work correctly.

set.seed(1234)
x2 <- rnorm(100, 1e-4, 2e-4)
for(i in 6:0)
    cat(i, try(mledist(x*10^i, "cauchy")$estimate, silent=TRUE), "\n")
        
    
}

\keyword{ distribution }% at least one, from doc/KEYWORDS
