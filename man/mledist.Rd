\name{mledist}
\alias{mledist}
\alias{mle}
\title{ Maximum likelihood fit of univariate distributions}

\description{
 Fit of univariate distributions using maximum likelihood for censored or non censored data.  
}

\usage{
mledist(data, distr, start=NULL, fix.arg=NULL, optim.method="default", 
    lower=-Inf, upper=Inf, custom.optim=NULL, ...)
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
    
    NB: if your data are particularly low or high, a scaling may be needed before the optimization process. See example (7).
}

\value{ 
    \code{mledist} returns a list with following components,
    \item{ estimate }{ the parameter estimates }
    \item{ convergence }{ an integer code for the convergence of \code{optim} defined as below
    or defined by the user in the user-supplied optimization function. 
    
        \code{0} indicates successful convergence.
        
        \code{1} indicates that the iteration limit of \code{optim} has been reached.

        \code{10} indicates degeneracy of the Nealder-Mead simplex.

        \code{100} indicates that \code{optim} encountered an internal error.
        }
        
    \item{ loglik }{ the log-likelihood }
    \item{ hessian }{ a symmetric matrix computed by \code{optim} as an estimate of the Hessian
        at the solution found or computed in the user-supplied optimization function.
        It is used in \code{fitdist} to estimate standard errors. }
    \item{ optim.function }{ the name of the optimization function used for maximum likelihood  }
}


\seealso{ 
    \code{\link{mmedist}}, \code{\link{qmedist}}, \code{\link{fitdist}},\code{\link{fitdistcens}}, \code{\link{optim}},
    \code{\link{bootdistcens}} and \code{\link{bootdist}}.
}

\references{ 
    Venables W.N. and Ripley B.D. (2002) Modern applied statistics with S. Springer, New York, pp. 435-446.
}


\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vetagro-sup.fr} and
Christophe Dutang
}

\examples{

# (1) basic fit of a normal distribution with maximum likelihood estimation
#

x1<-c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
mledist(x1,"norm")

# (2) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view dedicated to probability distributions

dgumbel<-function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
mledist(x1,"gumbel",start=list(a=10,b=5))

# (3) fit a discrete distribution (Poisson)
#

x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
mledist(x2,"pois")
mledist(x2,"nbinom")

# (4) fit a finite-support distribution (beta)
#

x3<-c(0.80,0.72,0.88,0.84,0.38,0.64,0.69,0.48,0.73,0.58,0.81,
0.83,0.71,0.75,0.59)
mledist(x3,"beta")


# (5) fit frequency distributions on USArrests dataset.
#

x4 <- USArrests$Assault
mledist(x4, "pois")
mledist(x4, "nbinom")

# (6) fit a continuous distribution (Gumbel) to censored data.
#

d1<-data.frame(
left=c(1.73,1.51,0.77,1.96,1.96,-1.4,-1.4,NA,-0.11,0.55,0.41,
    2.56,NA,-0.53,0.63,-1.4,-1.4,-1.4,NA,0.13),
right=c(1.73,1.51,0.77,1.96,1.96,0,-0.7,-1.4,-0.11,0.55,0.41,
    2.56,-1.4,-0.53,0.63,0,-0.7,NA,-1.4,0.13))
mledist(d1,"norm")

dgumbel<-function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel<-function(q,a,b) exp(-exp((a-q)/b))
mledist(d1,"gumbel",start=list(a=0,b=2),optim.method="Nelder-Mead")


# (7) scaling problem
#

x <- c(-0.00707717, -0.000947418, -0.00189753, 
-0.000474947, -0.00190205, -0.000476077, 0.00237812, 0.000949668, 
0.000474496, 0.00284226, -0.000473149, -0.000473373, 0, 0, 0.00283688, 
-0.0037843, -0.0047506, -0.00238379, -0.00286807, 0.000478583, 
0.000478354, -0.00143575, 0.00143575, 0.00238835, 0.0042847, 
0.00237248, -0.00142281, -0.00142484, 0, 0.00142484, 0.000948767, 
0.00378609, -0.000472478, 0.000472478, -0.0014181, 0, -0.000946522, 
-0.00284495, 0, 0.00331832, 0.00283554, 0.00141476, -0.00141476, 
-0.00188947, 0.00141743, -0.00236351, 0.00236351, 0.00235794, 
0.00235239, -0.000940292, -0.0014121, -0.00283019, 0.000472255, 
0.000472032, 0.000471809, -0.0014161, 0.0014161, -0.000943842, 
0.000472032, -0.000944287, -0.00094518, -0.00189304, -0.000473821, 
-0.000474046, 0.00331361, -0.000472701, -0.000946074, 0.00141878, 
-0.000945627, -0.00189394, -0.00189753, -0.0057143, -0.00143369, 
-0.00383326, 0.00143919, 0.000479272, -0.00191847, -0.000480192, 
0.000960154, 0.000479731, 0, 0.000479501, 0.000958313, -0.00383878, 
-0.00240674, 0.000963391, 0.000962464, -0.00192586, 0.000481812, 
-0.00241138, -0.00144963)


#only i == 0, no scaling, should not converge.
for(i in 6:0)
	cat(i, try(mledist(x*10^i, "cauchy")$estimate, silent=TRUE), "\n")
		
	
}

\keyword{ distribution }% at least one, from doc/KEYWORDS
