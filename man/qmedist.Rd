\name{qmedist}
\alias{qmedist}
\alias{qme}
\title{ Quantile matching fit of univariate distributions}

\description{
 Fit of univariate distribution by matching quantiles for non censored data.  
}

\usage{
qmedist(data, distr, probs, start=NULL, fix.arg=NULL, qtype=7, 
    optim.method="default", lower=-Inf, upper=Inf, custom.optim=NULL, ...)

}

\arguments{
\item{data}{ A numeric vector for non censored data. }
\item{distr}{ A character string \code{"name"} naming a distribution  for which the corresponding 
    quantile function 
    \code{qname} and the corresponding density distribution \code{dname} must be classically defined. }
\item{probs}{A numeric vector of the probabilities for which the quantile matching is done.
    The length of this vector must be equal
    to the number of parameters to estimate.}  
\item{start}{ A named list giving the initial values of parameters of the named distribution. 
    This argument may be omitted for some distributions for which reasonable 
    starting values are computed (see details).  }
\item{fix.arg}{ An optional named list giving the values of parameters of the named distribution
    that must kept fixed rather than estimated. }
\item{qtype}{The quantile type used by the R \code{\link{quantile}} function to 
    compute the empirical quantiles, (default 7 corresponds to the default quantile method in R).}
\item{optim.method}{ \code{"default"} or optimization method to pass to \code{\link{optim}}.  }
\item{lower}{ Left bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{upper}{ Right bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{custom.optim}{a function carrying the optimization.}
\item{\dots}{further arguments passed to the \code{optim} or \code{custom.optim} function.}
}

\details{
    The \code{qmedist} function carries out the quantile matching numerically, by minimization of the
    sum of squared differences between observed and theoretical quantiles.
    The optimization process is the same as \code{\link{mledist}}, see the 'details' section 
    of  \code{\link{mledist}}.

    This function is not intended to be called directly but is internally called in
    \code{\link{fitdist}} and \code{\link{bootdist}}.
}

\value{ 
    \code{qmedist} returns a list with following components,
    \item{ estimate }{ the parameter estimates.}
    \item{ convergence }{ an integer code for the convergence of \code{optim} defined as below
    or defined by the user in the user-supplied optimization function. 
    
        \code{0} indicates successful convergence.
        
        \code{1} indicates that the iteration limit of \code{optim} has been reached.

        \code{10} indicates degeneracy of the Nealder-Mead simplex.

        \code{100} indicates that \code{optim} encountered an internal error.
        }
    \item{ value }{the value of the statistic distance corresponding to estimate.}    
    \item{ hessian }{ a symmetric matrix computed by \code{optim} as an estimate of the Hessian
        at the solution found or computed in the user-supplied optimization function. }
    \item{ probs }{ the probability vector on which quantiles are matched. }   
    \item{ optim.function }{ the name of the optimization function used.  }
    \item{ loglik }{ the log-likelihood. }
}


\seealso{ 
    \code{\link{mmedist}}, \code{\link{mledist}}, \code{\link{fitdist}} for other estimation methods and 
    \code{\link{quantile}} for empirical quantile estimation in R.
}

\references{ 
    Klugman et al. (2004) Loss Models: From Data to Decisions, 2nd edition. 
    Wiley Series in Probability and Statistics, p. 331.
}


\author{ 
Christophe Dutang and Marie Laure Delignette-Muller.
}

\examples{

# (1) basic fit of a normal distribution 
#

x1<-c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
qmedist(x1, "norm", probs=c(1/3, 2/3))


# (2) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view dedicated to probability distributions

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
qgumbel <- function(p, a, b) a - b*log(-log(p))
qmedist(x1, "gumbel", probs=c(1/3, 2/3), start=list(a=10,b=5))

# (3) fit a discrete distribution (Poisson)
#

x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
qmedist(x2, "pois", probs=1/2)
qmedist(x2, "nbinom", probs=c(1/3, 2/3))

# (4) fit a finite-support distribution (beta)
#

x3<-c(0.80,0.72,0.88,0.84,0.38,0.64,0.69,0.48,0.73,0.58,0.81,
0.83,0.71,0.75,0.59)
qmedist(x3, "beta", probs=c(1/3, 2/3))


# (5) fit frequency distributions on USArrests dataset.
#

x4 <- USArrests$Assault
qmedist(x4, "pois", probs=1/2)
qmedist(x4, "nbinom", probs=c(1/3, 2/3))



}
\keyword{ distribution }% at least one, from doc/KEYWORDS
