\name{quantiledistcens}
\alias{quantiledistcens}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Estimation of quantiles from a distribution fitted on censored data}

\description{
 Uses the fit of a distribution to censored data to estimate quantiles for given probabilities, 
 using non parametric bootstrap to calculate a 95 percent confidence interval for each estimated quantile (optional).  
}

\usage{
quantiledistcens(f, probs = seq(0.1,0.9,0.1), bootstrap = TRUE, 
                bootstrap.arg = list(niter=1001))
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{f}{An object of class 'fitdistcens', result of the function \code{fitdistcens}.}
\item{probs}{numeric vector of probabilities with values in [0,1].}
\item{bootstrap}{If \code{TRUE}, bootstrap is performed to calculate quantile confidence intervals. }
\item{bootstrap.arg}{list of arguments to be passed to \code{\link{bootdistcens}}  }
}

\details{
    Quantiles of the parametric distribution are calculated for
    each probability specified in \code{probs}, using the estimated parameters given by the 'fitdist' result.
    If \code{bootstrap} is \code{TRUE}, a bootstrap procedure is then performed using \code{\link{bootdistcens}}
    on the 'fitdist' result, and quantiles of the parametric distribution are calculated for
    each probability specified in \code{probs}, using each set of boostrapped values of parameters. 
    The 95 percent confidence intervals of quantiles (2.5 and 97.5 percentiles) are calculated.
    
}

\value{ 
    \code{quantiledist} returns a list with 4 components,
    \item{ quantiles }{ a dataframe containing the estimated quantiles for each probability value specified in 
        the argument \code{probs} (one row, and as many columns as values in \code{probs}}
    \item{ resbootdistcens }{ the output of the function \code{\link{bootdistcens}} applied to \code{f} or NULL if 
        \code{bootstrap} is \code{FALSE}}
    \item{bootquant}{ a data frame containing the bootstraped values for each quantile
         (many rows, as specified in the call to \code{\link{bootdist}},
         and as many columns as values in \code{probs}, or NULL if 
        \code{bootstrap} is \code{FALSE} }
    \item{quantCI}{ bounds of the  95 percent confidence intervals for each quantile 
        (two rows, the first one for the 2.5 percentile and the second one for the 97.5 percentile,
         and as many columns as values in \code{probs}, or NULL if 
        \code{bootstrap} is \code{FALSE}}
}
        

\seealso{ 
    \code{\link{fitdistcens}}, \code{\link{bootdistcens}}.
}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} 
}

%\note{  }

\examples{
# (1) Fit of a normal distribution on censored data
# followed by estimation of quantiles of 
# the fitted distribution for various probabilities
#
d1<-data.frame(
left=c(1.73,1.51,0.77,1.96,1.96,-1.4,-1.4,NA,-0.11,0.55,0.41,
    2.56,NA,-0.53,0.63,-1.4,-1.4,-1.4,NA,0.13),
right=c(1.73,1.51,0.77,1.96,1.96,0,-0.7,-1.4,-0.11,0.55,0.41,
    2.56,-1.4,-0.53,0.63,0,-0.7,NA,-1.4,0.13))
f1n<-fitdistcens(d1, "norm")
f1n
q1 <- quantiledistcens(f1n, probs=c(0.05,0.1,0.2,0.5),bootstrap = FALSE)
q1b <- quantiledistcens(f1n, probs=c(0.05,0.1,0.2,0.5),bootstrap = TRUE)
q1b

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
