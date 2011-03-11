\name{mgedist}
\alias{mgedist}
\alias{mge}
\title{ Maximum goodness-of-fit fit of univariate continuous distributions}

\description{
 Fit of univariate continuous distribution by maximizing a goodness-of-fit distance for non censored data.  
}

\usage{
mgedist(data, distr, gof="CvM", start=NULL, fix.arg=NULL, 
    optim.method="default", lower=-Inf, upper=Inf, custom.optim=NULL, ...)

}

\arguments{
\item{data}{ A numeric vector for non censored data. }
\item{distr}{ A character string \code{"name"} naming a distribution  for which the corresponding 
    quantile function 
    \code{qname} and the corresponding density distribution \code{dname} must be classically defined. }
\item{gof}{A character string coding for the name of the goodness-of-fit distance used : 
    "CvM" for Cramer-von Mises distance,"KS" for Kolmogorov-Smirnov distance,
    "AD" for Anderson-Darling distance, "ADR", "ADL", "AD2R", "AD2L" and "AD2" for variants of
    Anderson-Darling distance described by Luceno (2006).}  
\item{start}{ A named list giving the initial values of parameters of the named distribution. 
    This argument may be omitted for some distributions for which reasonable 
    starting values are computed (see details).  }
\item{fix.arg}{ An optional named list giving the values of parameters of the named distribution
    that must kept fixed rather than estimated. }
\item{optim.method}{ \code{"default"} or optimization method to pass to \code{\link{optim}}.  }
\item{lower}{ Left bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{upper}{ Right bounds on the parameters for the \code{"L-BFGS-B"} method (see \code{\link{optim}}).  }
\item{custom.optim}{a function carrying the optimization.}
\item{\dots}{further arguments passed to the \code{optim} or \code{custom.optim} function.}
}

\details{
    The \code{mgedist} function numerically maximizes the goodness-of-fit distance coded by the argument 
    \code{gof}. One may use the classical Cramer-von Mises distance (\code{"CvM"}), the classical
    Kolmogorov-Smirnov distance (\code{"KS"}), the classical Anderson-Darling distance (\code{"AD"})
    which gives more weight to the tails of the distribution,
    or one of the variants of this last distance proposed by Luceno (2006). The right-tail AD (\code{"ADR"})
    gives more weight only to the right tail, the left-tail AD (\code{"ADL"})
    gives more weight only to the left tail. Either of the tails, or both of them, can receive even larger
    weights by using second order Anderson-Darling Statistics (using \code{"AD2R"}, \code{"AD2L"} or \code{"AD2"}).
    The minimization process is the same as \code{\link{mledist}}, see the 'details' section 
    of  \code{\link{mledist}}.

    This function is not intended to be called directly but is internally called in
    \code{\link{fitdist}} and \code{\link{bootdist}}.
    
    This function is intended to be used only with continuous distributions.
}

\value{ 
    \code{mgedist} returns a list with following components,
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
    \item{ gof }{ the code of the goodness-of-fit distance maximized. }   
    \item{ optim.function }{ the name of the optimization function used.  }
    \item{ loglik }{ the log-likelihood. }
}


\seealso{ 
    \code{\link{mmedist}}, \code{\link{mledist}}, \code{\link{qmedist}}, 
    \code{\link{fitdist}} for other estimation methods.
}

\references{ 
    Luceno, A. 2006. Fitting the generalized Pareto distribution to data using
    maximum goodness-of-fit estimators. \emph{Computational Statistics and Data Analysis}, \bold{51}, 904-917.
}


\author{ 
Marie Laure Delignette-Muller.
}

\examples{

# (1) Fit of a Weibull distribution to serving size data by maximum 
# goodness-of-fit estimation using all the distances available
# 

data(groundbeef)
serving <- groundbeef$serving
mgedist(serving,"weibull",gof="CvM")
mgedist(serving,"weibull",gof="KS")
mgedist(serving,"weibull",gof="AD")
mgedist(serving,"weibull",gof="ADR")
mgedist(serving,"weibull",gof="ADL")
mgedist(serving,"weibull",gof="AD2R")
mgedist(serving,"weibull",gof="AD2L")
mgedist(serving,"weibull",gof="AD2")

# (2) Fit of a uniform distribution using Cramer-von Mises or
# Kolmogorov-Smirnov distance
# 

u <- runif(100,min=5,max=10)
mgedist(u,"unif",gof="CvM")
mgedist(u,"unif",gof="KS")


}
\keyword{ distribution }% at least one, from doc/KEYWORDS
