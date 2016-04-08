\name{mledist}
\alias{mledist}
\alias{mle}
\title{ Maximum likelihood fit of univariate distributions}

\description{
 Fit of univariate distributions using maximum likelihood for censored or non censored data.  
}

\usage{
mledist(data, distr, start = NULL, fix.arg = NULL, optim.method = "default", 
    lower = -Inf, upper = Inf, custom.optim = NULL, weights = NULL, silent = TRUE, \dots)
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
     for which the corresponding density function \code{dname} 
     and the corresponding distribution function \code{pname} must be classically defined. }
\item{start}{A named list giving the initial values of parameters of the named distribution
    or a function of data computing initial values and returning a named list.
    This argument may be omitted (default) for some distributions for which reasonable 
    starting values are computed (see details). }
\item{fix.arg}{An optional named list giving the values of fixed parameters of the named distribution
    or a function of data computing (fixed) parameter values and returning a named list. 
    Parameters with fixed value are thus NOT estimated by this maximum likelihood procedure.}
\item{optim.method}{ \code{"default"} (see details) or an optimization method 
    to pass to \code{\link{optim}}.  }
\item{lower}{Left bounds on the parameters for the \code{"L-BFGS-B"} method 
    (see \code{\link{optim}}).  }
\item{upper}{ Right bounds on the parameters for the \code{"L-BFGS-B"} method 
    (see \code{\link{optim}}).  }
\item{custom.optim}{a function carrying the MLE optimisation (see details).}
\item{weights}{an optional vector of weights to be used in the fitting process. 
    Should be \code{NULL} or a numeric vector. If non-\code{NULL}, 
    weighted MLE is used, otherwise ordinary MLE.}
\item{silent}{A logical to remove or show warnings when bootstraping.}    
\item{\dots}{further arguments passed to the \code{optim} or \code{custom.optim} function.}
}

\details{
    This function is not intended to be called directly but is internally called in
    \code{\link{fitdist}} and \code{\link{bootdist}} when used with the maximum likelihood method
    and \code{\link{fitdistcens}} and \code{\link{bootdistcens}}. 

    It is assumed that the \code{distr} argument specifies the distribution by the
    probability density function and the cumulative distribution function (d, p).
    The quantile function and the random generator function (q, r) may be
    needed by other function such as \code{\link{mmedist}}, \code{\link{qmedist}}, \code{\link{mgedist}}, 
    \code{\link{fitdist}},\code{\link{fitdistcens}}, \code{\link{bootdistcens}} and \code{\link{bootdist}}.
  
    For the following named distributions, reasonable starting values will 
    be computed if \code{start} is omitted (i.e. \code{NULL}) : \code{"norm"}, \code{"lnorm"},
    \code{"exp"} and \code{"pois"}, \code{"cauchy"}, \code{"gamma"}, \code{"logis"},
    \code{"nbinom"} (parametrized by mu and size), \code{"geom"}, \code{"beta"}, \code{"weibull"}
    from the \code{stats} package; 
    \code{"invgamma"}, \code{"llogis"}, \code{"invweibull"}, \code{"pareto1"}, \code{"pareto"}
    from the \code{actuar} package. 
    Note that these starting values may not be good enough if the fit is poor. 
    The function uses a closed-form formula to fit the uniform distribution.
    If \code{start} is a list, then it should be a named list with the same names as in
    the d,p,q,r functions of the chosen distribution.
    If \code{start} is a function of data, then the function should return a named list with the same names as in
    the d,p,q,r functions of the chosen distribution.
    
    The \code{mledist} function allows user to set a fixed values for some parameters.
    As for \code{start}, if \code{fix.arg} is a list, then it should be a named list with the same names as in
    the d,p,q,r functions of the chosen distribution.
    If \code{fix.arg} is a function of data, then the function should return a named list with the 
    same names as in the d,p,q,r functions of the chosen distribution.

    When \code{custom.optim=NULL} (the default), maximum likelihood estimations 
    of the distribution parameters are computed with the R base \code{\link{optim}}.
    Direct optimization of the log-likelihood is performed (using \code{\link{optim}}) by default with 
    the \code{"Nelder-Mead"} method for distributions characterized by more than one parameter and
    the \code{"BFGS"} method for distributions characterized by only one parameter, or 
    with the method specified in the argument \code{"optim.method"} if not \code{"default"}.
    Box-constrainted optimization may be used with the method \code{"L-BFGS-B"}, using 
    the constraints on parameters specified in arguments \code{lower} and \code{upper}. 
    If non-trivial bounds are supplied, this method will be automatically selected, with a 
    warning.
    When errors are raised by \code{optim}, it's a good idea to start by adding traces during
    the optimization process by adding \code{control=list(trace=1, REPORT=1)}.
    
    If \code{custom.optim} is not \code{NULL}, then the user-supplied function is used 
    instead of the R base \code{optim}. The \code{custom.optim} must have (at least) 
    the following arguments
    \code{fn} for the function to be optimized, \code{par} for the initialized parameters.
    Internally the function to be optimized will also have other arguments, 
    such as \code{obs} with observations and \code{ddistname}
    with distribution name for non censored data (Beware of potential conflicts with optional
    arguments of \code{custom.optim}). It is assumed that \code{custom.optim} should carry 
    out a MINIMIZATION. 
    Finally, it should return at least the following components \code{par} for the estimate, 
    \code{convergence} for the convergence code, \code{value} for \code{fn(par)}, 
    \code{hessian} and \code{counts} for the number of calls (function and gradient), 
    see the returned value of \code{\link{optim}}. 
    See examples in \code{\link{fitdist}} and \code{\link{fitdistcens}}.
    
    Optionally, a vector of \code{weights} can be used in the fitting process. 
    By default (when \code{weigths=NULL}), ordinary MLE is carried out, otherwise 
    the specified weights are used to balance the log-likelihood contributions.
    It is not yet possible to take into account weighths in functions plotdist, plotdistcens,
    plot.fitdist, 
    plot.fitdistcens, cdfcomp, cdfcompcens, denscomp, ppcomp, qqcomp, gofstat and descdist 
    (developments planned in the future).

    
    NB: if your data values are particularly small or large, a scaling may be needed 
    before the optimization process. See example (7).
}

\value{ 
    \code{mledist} returns a list with following components,
    \item{estimate}{ the parameter estimates.}
    \item{convergence}{ an integer code for the convergence of \code{optim} defined as below
    or defined by the user in the user-supplied optimization function. 
    
        \code{0} indicates successful convergence.
        
        \code{1} indicates that the iteration limit of \code{optim} has been reached.

        \code{10} indicates degeneracy of the Nealder-Mead simplex.

        \code{100} indicates that \code{optim} encountered an internal error.
        }
        
    \item{loglik}{ the log-likelihood value.}
    \item{hessian}{ a symmetric matrix computed by \code{optim} as an estimate of the Hessian
        at the solution found or computed in the user-supplied optimization function.
        It is used in \code{fitdist} to estimate standard errors. }
    \item{optim.function}{ the name of the optimization function used for maximum likelihood.}
    \item{fix.arg}{ the named list giving the values of parameters of the named distribution
        that must kept fixed rather than estimated by maximum likelihood or \code{NULL} 
        if there are no such parameters. }
    \item{optim.method}{when \code{optim} is used, the name of the algorithm used, 
    \code{NULL} otherwise.}
    \item{fix.arg.fun}{the function used to set the value of \code{fix.arg} or \code{NULL}.}
    \item{weights}{the vector of weigths used in the estimation process or \code{NULL}.}
     \item{counts}{A two-element integer vector giving the number of calls 
    to the log-likelihood function and its gradient respectively. 
    This excludes those calls needed to compute the Hessian, if requested, 
    and any calls to log-likelihood function to compute a finite-difference
    approximation to the gradient. \code{counts} is returned by \code{optim}
    or the user-supplied function or set to \code{NULL}.}
    \item{optim.message}{A character string giving any additional information 
    returned by the optimizer, or \code{NULL}.}
}


\seealso{ 
    \code{\link{mmedist}}, \code{\link{qmedist}}, \code{\link{mgedist}}, 
    \code{\link{fitdist}},\code{\link{fitdistcens}}, \code{\link{optim}},
    \code{\link{bootdistcens}} and \code{\link{bootdist}}.
}

\references{ 

Venables WN and Ripley BD (2002), \emph{Modern applied statistics with S}.
Springer, New York, pp. 435-446.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.
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
# the simulated dataset (below) has particularly small values, 
# hence without scaling (10^0),
# the optimization raises an error. The for loop shows how scaling by 10^i
# for i=1,...,6 makes the fitting procedure work correctly.

set.seed(1234)
x2 <- rnorm(100, 1e-4, 2e-4)
for(i in 6:0)
    cat(i, try(mledist(x*10^i, "cauchy")$estimate, silent=TRUE), "\n")
        
 
# (17) small example for the zero-modified geometric distribution
#

dzmgeom <- function(x, p1, p2) p1 * (x == 0) + (1-p1)*dgeom(x-1, p2) #pdf
x2 <- c(2,  4,  0, 40,  4, 21,  0,  0,  0,  2,  5,  0,  0, 13,  2) #simulated dataset
initp1 <- function(x) list(p1=mean(x == 0)) #init as MLE
mledist(x2, "zmgeom", fix.arg=initp1, start=list(p2=1/2))

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
