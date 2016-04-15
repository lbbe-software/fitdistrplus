\name{prefit}
\alias{prefit}
\title{Pre-fitting procedure}

\description{
 Search good starting values
}

\usage{
prefit(data, distr, method = c("mle", "mme", "qme", "mge"), feasible.par,
  memp=NULL, order=NULL, probs=NULL, qtype=7, fix.arg=NULL, lower, 
  upper, weights=NULL, silent=TRUE, \dots)


}

\arguments{
\item{data}{A numeric vector.}
\item{distr}{A character string \code{"name"} naming a distribution for which the corresponding
    density function \code{dname}, the corresponding distribution function \code{pname} and the 
    corresponding quantile function \code{qname} must be defined, or directly the density function.}
\item{method}{A character string coding for the fitting method: 
    \code{"mle"} for 'maximum likelihood estimation', \code{"mme"} for 'moment matching estimation', 
    \code{"qme"} for 'quantile matching estimation' and \code{"mge"} for 'maximum goodness-of-fit estimation'.}
\item{feasible.par}{A named list giving the initial values of parameters of the named distribution
    or a function of data computing initial values and returning a named list.
    This argument may be omitted (default) for some distributions for which reasonable 
    starting values are computed (see the 'details' section of  \code{\link{mledist}}). 
    It may not be into account for closed-form formulas.}
\item{order}{A numeric vector for the moment order(s). The length of this vector must be equal
to the number of parameters to estimate.}
\item{memp}{A function implementing empirical moments, raw or centered but has to be consistent with 
    \code{distr} argument (and \code{weights} argument).}
\item{probs}{A numeric vector of the probabilities for which the quantile matching is done.
    The length of this vector must be equal to the number of parameters to estimate.}  
\item{qtype}{The quantile type used by the R \code{\link{quantile}} function to 
    compute the empirical quantiles, (default 7 corresponds to the default quantile method in R).}
    
\item{fix.arg}{An optional named list giving the values of fixed parameters of the named distribution
    or a function of data computing (fixed) parameter values and returning a named list. 
    Parameters with fixed value are thus NOT estimated by this maximum likelihood procedure.
    The use of this argument is not possible if \code{method="mme"} and a closed-form formula is used.}
\item{weights}{an optional vector of weights to be used in the fitting process. 
    Should be \code{NULL} or a numeric vector. If non-\code{NULL}, 
    weighted MLE is used, otherwise ordinary MLE.}
\item{silent}{A logical to remove or show warnings.}
\item{lower}{Lower bounds on the parameters.}
\item{upper}{Upper bounds on the parameters.}
    
\item{\dots}{Further arguments to be passed to generic functions, or to one of the functions
    \code{"mledist"},  \code{"mmedist"}, \code{"qmedist"} or \code{"mgedist"}
    depending of the chosen method. See \code{\link{mledist}}, \code{\link{mmedist}}, \code{\link{qmedist}}, 
    \code{\link{mgedist}} for details on parameter estimation.}
}

\details{

Search good starting values   
    
    
}

\value{ 
    A named list.
}

\seealso{
    See \code{\link{mledist}}, \code{\link{mmedist}}, \code{\link{qmedist}}, 
    \code{\link{mgedist}} for details on parameter estimation.
    See \code{\link{fitdist}} for the main procedure.
}

\references{ 
Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} 
and Christophe Dutang.
}

\examples{

# (1) fit of a gamma distribution by maximum likelihood estimation
#
x <- rgamma(1e3, 5/2, 7/2)

prefit(x, "gamma", "mle", list(shape=3, scale=3), lower=-Inf, upper=Inf)

}
\keyword{ distribution }
