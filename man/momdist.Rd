\name{momdist}
\alias{momdist}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fitting of univariate distributions by matching moments}

\description{
 Fits a univariate distribution by matching moments.  
}

\usage{
momdist(data, distr)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{data}{ A numeric vector.  }
\item{distr}{ A character string \code{"name"} naming a distribution or directly the density function
                dname.
                The estimated values of the distribution parameters are provided only for the following
                distributions : \code{"norm"}, \code{"lnorm"},
                \code{"exp"} and \code{"pois"}, \code{"gamma"}, \code{"logis"},
                \code{"nbinom"} , \code{"geom"}, \code{"beta"} and \code{"unif"}.
            }
}

\details{ For distributions characterized by one parameter (\code{"geom"}, \code{"pois"} and \code{"exp"}),
         this parameter is simply 
        estimated by matching theoretical and observed means, and for distributions characterized by 
        two parameters, these parameters are estimated by matching theoretical and observed means
        and variances (Vose, 2000). 
        
        This function is not intended to be called directly but is internally called in
        \code{\link{fitdist}} and \code{\link{bootdist}} when used with the matching moments method. 
     
}

\value{ 
    \code{momdist} returns the named parameter or a named vector of parameters.
}


\seealso{ 
    \code{\link{mledist}}, \code{\link{fitdist}} and \code{\link{bootdist}}.
}

\references{ 
    Vose D (2000) Risk analysis, a quantitative guide. John Wiley & Sons Ltd, Chischester, England, pp. 99-143.
    Evans M, Hastings N and Peacock B (2000) Statistical distributions. John Wiley and Sons Inc.
}


\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vet-lyon.fr} 
}

%\note{  }

\examples{
x1<-c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
momdist(x1,"norm")

x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
momdist(x2,"pois")

x3<-c(0.80,0.72,0.88,0.84,0.38,0.64,0.69,0.48,0.73,0.58,0.81,
0.83,0.71,0.75,0.59)
momdist(x3,"beta")
}

\keyword{ distribution }% at least one, from doc/KEYWORDS
