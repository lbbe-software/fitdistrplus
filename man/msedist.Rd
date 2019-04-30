\name{msedist}
\alias{msedist}
\alias{mse}
\title{Maximum spacing estimation of univariate distributions}

\description{
 Fit of univariate distribution by maximizing (log) spacings for non censored data.  
}

\usage{
msedist(data, distr, phidiv="KL", power.phidiv=NULL, start = NULL, fix.arg = NULL, 
  optim.method = "default", lower = -Inf, upper = Inf, custom.optim = NULL, 
  weights=NULL, silent = TRUE, gradient = NULL, checkstartfix=FALSE, \dots)

}

\arguments{
\item{data}{ A numeric vector for non censored data. }
\item{distr}{ A character string \code{"name"} naming a distribution  for which the corresponding 
    quantile function 
    \code{qname} and the corresponding density distribution \code{dname} must be classically defined. }
\item{phidiv}{A character string coding for the name of the phi-divergence used : 
    \code{"KL"} for Kullback-Leibler information (corresponds to classic maximum spacing estimation),
    \code{"J"} for Jeffreys' divergence, \code{"R"} for Renyi's divergence,
    \code{"H"} for Hellinger distance, \code{"V"} for Vajda's measure of information, see details.}  
\item{power.phidiv}{If relevant, a numeric for the power used in some phi-divergence : 
    should be \code{NULL} when \code{phidiv="KL"} or \code{phidiv="J"} ,
    should be positive and different from 1 when \code{phidiv="R"},
    should be greater or equal to 1 when \code{phidiv="H"} or \code{phidiv="V"}, see details.}    
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
\item{custom.optim}{a function carrying the optimization.}
\item{weights}{an optional vector of weights to be used in the fitting process. 
    Should be \code{NULL} or a numeric vector with strictly positive integers 
    (typically the number of occurences of each observation). If non-\code{NULL}, 
    weighted MLE is used, otherwise ordinary MLE.}
\item{silent}{A logical to remove or show warnings when bootstraping.}    
\item{gradient}{A function to return the gradient of the gof distance for the \code{"BFGS"}, \code{"CG"} 
    and \code{"L-BFGS-B"} methods. If it is \code{NULL}, a finite-difference approximation will be used.}
\item{checkstartfix}{A logical to test starting and fixed values. Do not change it.}      
\item{\dots}{further arguments passed to the \code{\link{optim}}, \code{\link{constrOptim}} or \code{custom.optim} function.}
}

\details{
    The \code{msedist} function numerically maximizes a phi-divergence function of spacings,
    where spacings are the differences of the cumulative distribution function evaluated at
    the sorted dataset.
    The classical maximum spacing estimation (MSE) was introduced by Cheng and Amin (1986) 
    and Ranneby (1984) independently where the phi-diverence is the logarithm,
    see Anatolyev and Kosenok (2005) for a link between MSE and maximum likelihood estimation.
    
    MSE was generalized by Ranneby and Ekstrom (1997) by allowing different phi-divergence
    function. Generalized MSE maximizes
    \deqn{
      S_n(\theta)=\frac{1}{n+1}\sum_{i=1}^{n+1} \phi\left(F(x_{(i)}; \theta)-F(x_{(i-1)}; \theta) \right),
    }{
      S_n(\theta)= sum(i=1,..,n+1;  \phi(F(x_(i); \theta)-F(x_(i-1); \theta))) / (n+1),
    }
    where \eqn{F(;\theta)} is the parametric distribution function to be fitted,
    \eqn{\phi} is the phi-divergence function,
    \eqn{x_{(1)}<\dots<x_{(n)}}{x_(1)<\dots<x_(n)} is the sorted sample, 
    \eqn{x_{(0)}=-\infty}{x_(0)=-\infty}
    and  \eqn{x_{(n+1)}=+\infty}{x_(n+1)=+\infty}.
    The possible phi-divergence function is
    \itemize{
      \item Kullback-Leibler information (when \code{phidiv="KL"} and corresponds to classical MSE) 
      \deqn{\phi(x)=\log(x)}{\phi(x)=\log(x)}
      \item Jeffreys' divergence (when \code{phidiv="J"}) 
      \deqn{\phi(x)=(1-x)\log(x)}{\phi(x)=(1-x)\log(x)}
      \item Renyi's divergence (when \code{phidiv="R"} and \code{power.phidiv=alpha}) 
      \deqn{\phi(x)=x^\alpha\times\textrm{sign}(1-\alpha) \textrm{ with } \alpha>0, \alpha\neq 1
      }{
      \phi(x)=x^\alpha*sign(1-\alpha) with 0<\alpha<1 or 1 <\alpha
      } 
      \item Hellinger distance (when \code{phidiv="H"} and \code{power.phidiv=p})
      \deqn{\phi(x)=-|1-x^{1/p}|^p \textrm{ with } p\ge 1
      }{
      \phi(x)=-|1-x^(1/p)|^p with p\ge 1
      } 
      \item Vajda's measure of information (when \code{phidiv="V"} and \code{power.phidiv=beta})
      \deqn{\phi(x)=-|1-x|^\beta \textrm{ with } \beta\ge 1
      }{
      \phi(x)=-|1-x|^\beta with \beta\ge 1
      } 
    }
    
    
    
    The optimization process is the same as \code{\link{mledist}}, see the 'details' section 
    of that function.

    This function is not intended to be called directly but is internally called in
    \code{\link{fitdist}} and \code{\link{bootdist}}.
    
    This function is intended to be used only with non-censored data.
    
    NB: if your data values are particularly small or large, a scaling may be needed 
    before the optimization process, see \code{\link{mledist}}'s examples. 
}

\value{ 
    \code{msedist} returns a list with following components,
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
    \item{loglik}{ the log-likelihood value.}
    \item{phidiv}{The character string coding for the name of the phi-divergence used
    either \code{"KL"}, \code{"J"}, \code{"R"}, \code{"H"} or \code{"V"}.}  
    \item{power.phidiv}{Either \code{NULL} or a numeric for the power used in the phi-divergence.}    
}


\seealso{ 
    \code{\link{mmedist}}, \code{\link{mledist}}, \code{\link{qmedist}}, \code{\link{mgedist}}, 
    \code{\link{fitdist}} for other estimation methods.
}

\references{

Anatolyev, S., and Kosenok, G. (2005). \emph{An alternative to maximum likelihood
based on spacings}. Econometric Theory, 21(2), 472-476. 

Cheng, R.C.H. and N.A.K. Amin (1983) \emph{Estimating parameters in continuous univariate 
distributions with a shifted origin}. Journal of the Royal Statistical Society Series B 45, 394-403.

Ranneby, B. (1984) \emph{The maximum spacing method: An estimation method related to 
the maximum likelihood method}. Scandinavian Journal of Statistics 11, 93-112.

Ranneby, B. and Ekstroem, M. (1997). \emph{Maximum spacing estimates based
on different metrics}. Umea universitet.

}


\author{ 
Marie-Laure Delignette-Muller and Christophe Dutang.
}

\examples{

# (1) Fit of a Weibull distribution to serving size data by maximum 
# spacing estimation
# 

data(groundbeef)
serving <- groundbeef$serving
msedist(serving, "weibull")



# (2) Fit of an exponential distribution 
# 

set.seed(123)
x1 <- rexp(1e3)
#the convergence is quick
msedist(x1, "exp", control=list(trace=0, REPORT=1))

}
\keyword{ distribution }% at least one, from doc/KEYWORDS
