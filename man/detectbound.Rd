\name{detectbound}
\alias{detectbound}
\title{ Detect bounds for density function}

\description{
 Manual detection of bounds of parameter of a density function/  
}

\usage{
detectbound(distname, vstart, obs, fix.arg=NULL, echo=FALSE)
}

\arguments{
\item{distname}{ A character string \code{"name"} naming a distribution 
     for which the corresponding density function \code{dname} 
     must be classically defined. }
\item{vstart}{A named vector giving the initial values of parameters 
    of the named distribution. }
\item{obs}{ A numeric vector for non censored data. }
\item{fix.arg}{An optional named vector giving the values of fixed parameters of the named distribution. Default to \code{NULL}.}
\item{echo}{A logical to show some traces.}    
}

\details{
    This function manually tests the following bounds : -1, 0, and 1.
}

\value{ 
    \code{detectbound} returns a 2-row matrix with the lower bounds in the first
    row and the upper bounds in the second row.
}


\seealso{ 
    \code{\link{fitdist}}.
}

\references{ 


Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.
}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} and
Christophe Dutang.
}

\examples{

#case where the density returns a Not-an-Numeric value.
detectbound("exp", c(rate=3), 1:10)
detectbound("binom", c(size=3, prob=1/2), 1:10)
detectbound("nbinom", c(size=3, prob=1/2), 1:10)


}

\keyword{ distribution }% at least one, from doc/KEYWORDS
