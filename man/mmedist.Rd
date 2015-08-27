\name{mmedist}
\alias{mmedist}
\alias{mme}
\title{ Matching moment fit of univariate distributions}

\description{
 Fit of univariate distributions by matching moments (raw or centered) for non censored data.
}

\usage{
mmedist(data, distr, order, memp, start = NULL, fix.arg = NULL, optim.method = "default", 
  lower = -Inf, upper = Inf, custom.optim = NULL, weights = NULL, silent = TRUE, \dots) 

}
\arguments{
\item{data}{A numeric vector for non censored data.}
\item{distr}{A character string \code{"name"} naming a distribution (see 'details').}
\item{order}{A numeric vector for the moment order(s). The length of this vector must be equal
to the number of parameters to estimate.}
\item{memp}{A function implementing empirical moments, raw or centered but has to be consistent with 
    \code{distr} argument (and \code{weights} argument). 
    See details below.
    }
\item{start}{A named list giving the initial values of parameters of the named distribution
    or a function of data computing initial values and returning a named list.
    This argument may be omitted (default) for some distributions for which reasonable 
    starting values are computed (see the 'details' section of  \code{\link{mledist}}). }
\item{fix.arg}{An optional named list giving the values of fixed parameters of the named distribution
    or a function of data computing (fixed) parameter values and returning a named list. 
    Parameters with fixed value are thus NOT estimated.}
\item{optim.method}{ \code{"default"} or optimization method to pass to \code{\link{optim}}.  }
\item{lower}{ Left bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{upper}{ Right bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{custom.optim}{a function carrying the optimization .}
\item{weights}{an optional vector of weights to be used in the fitting process. 
    Should be \code{NULL} or a numeric vector. If non-\code{NULL}, 
    weighted MME is used, otherwise ordinary MME.}
\item{silent}{A logical to remove or show warnings when bootstraping.}        
\item{\dots}{further arguments passed to the \code{optim} or \code{custom.optim} function.} 
}

\details{ 
    The argument \code{distr} can be one of the base R distributions: \code{"norm"}, \code{"lnorm"},
    \code{"exp"} and \code{"pois"}, \code{"gamma"}, \code{"logis"},
    \code{"nbinom"} , \code{"geom"}, \code{"beta"} and \code{"unif"}.
    In that case, no other arguments than \code{data} and \code{distr} are
    required, because the estimate is computed by a closed-form formula.
    For distributions characterized by one parameter (\code{"geom"}, \code{"pois"} and \code{"exp"}),
    this parameter is simply estimated by matching theoretical and observed means, and for distributions 
    characterized by two parameters, these parameters are estimated by matching theoretical and observed 
    means and variances (Vose, 2000).
    Note that for these closed-form formula, \code{fix.arg} cannot be used and \code{start} is ignored.
    
    The argument \code{distr} can also be the distribution name
    as long as a corresponding \code{mdistr} function exists, e.g. \code{"pareto"} if \code{"mpareto"} exists.
    In that case arguments arguments \code{order} and \code{memp} 
    have to be supplied in order to carry out the matching numerically, by minimization of the
    sum of squared differences between observed and theoretical moments. 
    Optionnally other arguments can be supplied to control optimization (see the 'details' section of 
    \code{\link{mledist}} for details about arguments for the control of optimization). 
    In that case, \code{fix.arg} can be used and \code{start} is taken into account.
    
    For non closed-form estimators, \code{memp} must be provided to compute empirical moments.
    When \code{weights=NULL}, this function must have two arguments \code{x, order}: 
    \code{x} the numeric vector of the data and \code{order} the order of the moment.
    When \code{weights!=NULL}, this function must have three arguments \code{x, order, weights}: 
    \code{x} the numeric vector of the data, \code{order} the order of the moment,
    \code{weights} the numeric vector of weights. See examples below.

    Optionally, a vector of \code{weights} can be used in the fitting process. 
    By default (when \code{weigths=NULL}), ordinary MME is carried out, otherwise 
    the specified weights are used to compute (raw or centered) weighted moments.
    For closed-form estimators, weighted mean and variance are computed by 
    \code{\link[Hmisc:wtd.stats]{wtd.mean}} and \code{\link[Hmisc:wtd.stats]{wtd.var}}
    from the \code{Hmisc} package. When a numerical minimization
    is used, weighted are expected to be computed by the \code{memp} function.
    
    This function is not intended to be called directly but is internally called in
    \code{\link{fitdist}} and \code{\link{bootdist}} when used with the matching moments method. 
     
}

\value{ 
    \code{mmedist} returns a list with following components,
    \item{ estimate }{ the parameter estimates. }
    \item{ convergence }{ (if appropriate) an integer code for the convergence of \code{optim} defined as below
    or defined by the user in the user-supplied optimization function. 
    
        \code{0} indicates successful convergence.
        
        \code{1} indicates that the iteration limit of \code{optim} has been reached.

        \code{10} indicates degeneracy of the Nealder-Mead simplex.

        \code{100} indicates that \code{optim} encountered an internal error.
        }
    \item{loglik}{ the log-likelihood value.}
    \item{hessian}{ (if appropriate) a symmetric matrix computed by \code{optim} as an estimate of the Hessian
        at the solution found or computed in the user-supplied optimization function. }
    \item{optim.function}{ (if appropriate) the name of the optimization function. }
    \item{memp}{ (if appropriate) the empirical moment function. }
    \item{order}{ the order of the moment(s) matched.}        
    \item{method}{either \code{"closed formula"} or the name of the optimization method.}
    \item{fix.arg}{ the named list giving the values of parameters of the named distribution
        that must kept fixed rather than estimated by maximum likelihood or \code{NULL} 
        if there are no such parameters. }
    \item{optim.method}{when \code{optim} is used, the name of the algorithm used, 
    \code{NULL} otherwise.}        
    \item{weights}{the vector of weigths used in the estimation process or \code{NULL}.}
}


\seealso{ 
    \code{\link{mmedist}}, \code{\link{qmedist}}, \code{\link{fitdist}},\code{\link{fitdistcens}}, 
    \code{\link{optim}}, \code{\link{bootdistcens}} and \code{\link{bootdist}}.
}

\references{ 
    
Evans M, Hastings N and Peacock B (2000), \emph{Statistical distributions}. John Wiley and Sons Inc.

Vose D (2000), \emph{Risk analysis, a quantitative guide}.
John Wiley & Sons Ltd, Chischester, England, pp. 99-143.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} and 
Christophe Dutang.
}

\examples{

# (1) basic fit of a normal distribution with moment matching estimation
#

set.seed(1234)
n <- 100
x1 <- rnorm(n=n)
mmedist(x1, "norm")

#weighted
w <- c(rep(1, n/2), rep(10, n/2))
mmedist(x1, "norm", weights=w)$estimate


# (2) fit a discrete distribution (Poisson)
#

set.seed(1234)
x2 <- rpois(n=30,lambda = 2)
mmedist(x2, "pois")

# (3) fit a finite-support distribution (beta)
#

set.seed(1234)
x3 <- rbeta(n=100,shape1=5, shape2=10)
mmedist(x3, "beta")


# (4) fit a Pareto distribution
#

\dontrun{
  require(actuar)
  #simulate a sample
  x4  <-  rpareto(1000, 6, 2)

  #empirical raw moment
  memp  <-  function(x, order) mean(x^order)
  memp2 <- function(x, order, weights) sum(x^order * weights)/sum(weights)

  #fit by MME
  mmedist(x4, "pareto", order=c(1, 2), memp=memp, 
    start=list(shape=10, scale=10), lower=1, upper=Inf)
  #fit by weighted MME
  w <- rep(1, length(x4))
  w[x4 < 1] <- 2
  mmedist(x4, "pareto", order=c(1, 2), memp=memp2, weights=w,
    start=list(shape=10, scale=10), lower=1, upper=Inf)
}

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
