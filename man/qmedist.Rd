\name{qmedist}
\alias{qmedist}
\alias{qme}
\title{ Quantile matching fit of univariate distributions}

\description{
 Fit of univariate distribution by matching quantiles for non censored data.  
}

\usage{
qmedist(data, distr, probs, start = NULL, fix.arg = NULL, qtype = 7, 
    optim.method = "default", lower = -Inf, upper = Inf, 
    custom.optim = NULL, weights = NULL, silent = TRUE, gradient = NULL, 
    checkstartfix=FALSE, calcvcov=FALSE, \dots)

}

\arguments{
\item{data}{ A numeric vector for non censored data. }
\item{distr}{ A character string \code{"name"} naming a distribution  for which the corresponding 
    quantile function 
    \code{qname} and the corresponding density distribution \code{dname} must be classically defined. }
\item{probs}{A numeric vector of the probabilities for which the quantile matching is done.
    The length of this vector must be equal
    to the number of parameters to estimate.}  
\item{start}{A named list giving the initial values of parameters of the named distribution
    or a function of data computing initial values and returning a named list.
    This argument may be omitted (default) for some distributions for which reasonable 
    starting values are computed (see the 'details' section of  \code{\link{mledist}}). }
\item{fix.arg}{An optional named list giving the values of fixed parameters of the named distribution
    or a function of data computing (fixed) parameter values and returning a named list. 
    Parameters with fixed value are thus NOT estimated.}
\item{qtype}{The quantile type used by the R \code{\link{quantile}} function to 
    compute the empirical quantiles, (default 7 corresponds to the default quantile method in R).}
\item{optim.method}{ \code{"default"} or optimization method to pass to \code{\link{optim}}.  }
\item{lower}{ Left bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{upper}{ Right bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{custom.optim}{a function carrying the optimization.}
\item{weights}{an optional vector of weights to be used in the fitting process. 
    Should be \code{NULL} or a numeric vector with strictly positive integers 
    (typically the number of occurences of each observation). If non-\code{NULL}, 
    weighted QME is used, otherwise ordinary QME.}
\item{silent}{A logical to remove or show warnings when bootstraping.}        
\item{gradient}{A function to return the gradient of the squared difference for the \code{"BFGS"}, \code{"CG"} 
    and \code{"L-BFGS-B"} methods. If it is \code{NULL}, a finite-difference approximation will be used,
    see details.}
\item{checkstartfix}{A logical to test starting and fixed values. Do not change it.}     
\item{calcvcov}{A logical indicating if (asymptotic) covariance matrix is required.
  (currently ignored)}
\item{\dots}{further arguments passed to the \code{\link{optim}}, 
    \code{\link{constrOptim}} or \code{custom.optim} function.}

}

\details{
    The \code{qmedist} function carries out the quantile matching numerically, by minimization of the
    sum of squared differences between observed and theoretical quantiles.
    Note that for discrete distribution, the sum of squared differences is a step function and
    consequently, the optimum is not unique, see the FAQ.
    
    The optimization process is the same as \code{\link{mledist}}, see the 'details' section 
    of that function.
    
    Optionally, a vector of \code{weights} can be used in the fitting process. 
    By default (when \code{weigths=NULL}), ordinary QME is carried out, otherwise 
    the specified weights are used to compute weighted quantiles used in the squared differences.
    Weigthed quantiles are computed by \code{\link[Hmisc:wtd.stats]{wtdquantile}} from the \code{Hmisc} package.
    It is not yet possible to take into account weighths in functions \code{plotdist}, 
    \code{plotdistcens}, \code{plot.fitdist}, \code{plot.fitdistcens}, \code{cdfcomp}, 
    \code{cdfcompcens}, \code{denscomp}, \code{ppcomp}, \code{qqcomp}, \code{gofstat} 
    and \code{descdist}
    (developments planned in the future).


    This function is not intended to be called directly but is internally called in
    \code{\link{fitdist}} and \code{\link{bootdist}}.
}

\value{ 
    \code{qmedist} returns a list with following components,
    
    \item{estimate}{ the parameter estimates.}
    \item{convergence}{ an integer code for the convergence of \code{\link{optim}} defined as below
    or defined by the user in the user-supplied optimization function. 
        \code{0} indicates successful convergence.
        \code{1} indicates that the iteration limit of \code{\link{optim}} has been reached.
        \code{10} indicates degeneracy of the Nealder-Mead simplex.
        \code{100} indicates that \code{\link{optim}} encountered an internal error.
        }
    \item{value}{the minimal value reached for the criterion to minimize.}    
    \item{hessian}{ a symmetric matrix computed by \code{\link{optim}} as an estimate of the Hessian
        at the solution found or computed in the user-supplied optimization function. }
    \item{optim.function}{the name of the optimization function used for maximum likelihood.}
    \item{optim.method}{when \code{\link{optim}} is used, the name of the
    algorithm used, the field \code{method} of the \code{custom.optim} function
    otherwise.}
    \item{fix.arg}{the named list giving the values of parameters of the named distribution
        that must kept fixed rather than estimated by maximum likelihood or \code{NULL} 
        if there are no such parameters. }
    \item{fix.arg.fun}{the function used to set the value of \code{fix.arg} or \code{NULL}.}
    \item{weights}{the vector of weigths used in the estimation process or \code{NULL}.}
    \item{counts}{A two-element integer vector giving the number of calls 
    to the log-likelihood function and its gradient respectively. 
    This excludes those calls needed to compute the Hessian, if requested, 
    and any calls to log-likelihood function to compute a finite-difference
    approximation to the gradient. \code{counts} is returned by \code{\link{optim}}
    or the user-supplied function or set to \code{NULL}.}
    \item{optim.message}{A character string giving any additional information 
    returned by the optimizer, or \code{NULL}. To understand exactly the message,
    see the source code.}
    \item{loglik}{ the log-likelihood value. }
    \item{probs}{ the probability vector on which quantiles are matched. }   
}


\seealso{ 
    See \code{\link{mmedist}}, \code{\link{mledist}}, \code{\link{mgedist}},
    \code{\link{fitdist}} for other estimation methods and 
    \code{\link{quantile}} for empirical quantile estimation in R.
    
    Please visit the \href{https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html}{Frequently Asked Questions}.
}

\references{
Klugman SA, Panjer HH and Willmot GE (2012),
\emph{Loss Models: From Data to Decissions}, 4th edition.
Wiley Series in Statistics for Finance, Business and Economics, p. 253,
\doi{10.1198/tech.2006.s409}.


Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34, \doi{https://doi.org/10.18637/jss.v064.i04}.


}


\author{ 
Christophe Dutang and Marie Laure Delignette-Muller.

}

\examples{

# (1) basic fit of a normal distribution 
#

set.seed(1234)
x1 <- rnorm(n=100)
qmedist(x1, "norm", probs=c(1/3, 2/3))


# (2) defining your own distribution functions, here for the Gumbel 
# distribution for other distributions, see the CRAN task view dedicated 
# to probability distributions

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
qgumbel <- function(p, a, b) a - b*log(-log(p))
qmedist(x1, "gumbel", probs=c(1/3, 2/3), start=list(a=10,b=5))

# (3) fit a discrete distribution (Poisson)
#

set.seed(1234)
x2 <- rpois(n=30,lambda = 2)
qmedist(x2, "pois", probs=1/2)

# (4) fit a finite-support distribution (beta)
#

set.seed(1234)
x3 <- rbeta(n=100,shape1=5, shape2=10)
qmedist(x3, "beta", probs=c(1/3, 2/3))

# (5) fit frequency distributions on USArrests dataset.
#

x4 <- USArrests$Assault
qmedist(x4, "pois", probs=1/2)
qmedist(x4, "nbinom", probs=c(1/3, 2/3))

}
\keyword{ distribution }
