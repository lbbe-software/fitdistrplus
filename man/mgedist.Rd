\name{mgedist}
\alias{mgedist}
\alias{mge}
\title{ Maximum goodness-of-fit fit of univariate continuous distributions}

\description{
 Fit of univariate continuous distribution by maximizing goodness-of-fit (or minimizing distance) for non censored data.  
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
    The \code{mgedist} function numerically maximizes goodness-of-fit, 
    or minimizes a goodness-of-fit distance coded by the argument 
    \code{gof}. One may use one of the classical distances defined in Stephens (1986), 
    the Cramer-von Mises distance (\code{"CvM"}), the
    Kolmogorov-Smirnov distance (\code{"KS"}) or the Anderson-Darling distance (\code{"AD"})
    which gives more weight to the tails of the distribution,
    or one of the variants of this last distance proposed by Luceno (2006). The right-tail AD (\code{"ADR"})
    gives more weight only to the right tail, the left-tail AD (\code{"ADL"})
    gives more weight only to the left tail. Either of the tails, or both of them, can receive even larger
    weights by using second order Anderson-Darling Statistics (using \code{"AD2R"}, \code{"AD2L"} or \code{"AD2"}).
    The optimization process is the same as \code{\link{mledist}}, see the 'details' section 
    of  \code{\link{mledist}}.

    This function is not intended to be called directly but is internally called in
    \code{\link{fitdist}} and \code{\link{bootdist}}.
    
    This function is intended to be used only with continuous distributions.
    
    NB: if your data are particularly low or high, a scaling may be needed before the optimization process. See example (3).    
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

    Stephens MA (1986) Tests based on edf statistics. In Goodness-of-fit techniques (D'Agostino RB and
    Stephens MA, eds), Marcel dekker, New York, pp. 97-194.
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


# (3) scaling problem
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
	cat(i, try(mgedist(x*10^i,"cauchy")$estimate, silent=TRUE), "\n")
		
	

}
\keyword{ distribution }% at least one, from doc/KEYWORDS
