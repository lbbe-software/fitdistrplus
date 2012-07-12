\name{quantiledist}
\alias{quantiledist}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Estimation of quantiles from a distribution fitted on non-censored data}

\description{
 Uses the fit of a distribution to non censored data to estimate quantiles for given probabilities, 
 using parametric or nonparametric bootstrap to calculate a 95 percent confidence interval for 
 each estimated quantile (optional).  
}

\usage{
quantiledist(f, probs = seq(0.1,0.9,0.1), bootstrap = TRUE, 
                bootstrap.arg = list(bootmethod="param", niter=1001))
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{f}{An object of class 'fitdist', result of the function \code{fitdist}.}
\item{probs}{numeric vector of probabilities with values in [0,1].}
\item{bootstrap}{If \code{TRUE}, bootstrap is performed to calculate quantile confidence intervals. }
\item{bootstrap.arg}{list of arguments to be passed to \code{\link{bootdist}}  }
}

\details{
    Quantiles of the parametric distribution are calculated for
    each probability specified in \code{probs}, using the estimated parameters given by the 'fitdist' result.
    If \code{bootstrap} is \code{TRUE}, a bootstrap procedure is then performed using \code{\link{bootdist}}
    on the 'fitdist' result, and quantiles of the parametric distribution are calculated for
    each probability specified in \code{probs}, using each set of boostrapped values of parameters. 
    The 95 percent confidence intervals of quantiles (2.5 and 97.5 percentiles) are calculated.
    
}

\value{ 
    \code{quantiledist} returns a list with 4 components,
    \item{ quantiles }{ a dataframe containing the estimated quantiles for each probability value specified in 
        the argument \code{probs} (one row, and as many columns as values in \code{probs}}
    \item{ resbootdist }{ the output of the function \code{\link{bootdist}} applied to \code{f} or NULL if 
        \code{bootstrap} is \code{FALSE}}
    \item{bootquant }{ a data frame containing the bootstraped values for each quantile
         (many rows, as specified in the call to \code{\link{bootdist}},
         and as many columns as values in \code{probs}, or NULL if 
        \code{bootstrap} is \code{FALSE} }
    \item{quantCI}{ bounds of the  95 percent confidence intervals for each quantile 
        (two rows, the first one for the 2.5 percentile and the second one for the 97.5 percentile,
         and as many columns as values in \code{probs}, or NULL if 
        \code{bootstrap} is \code{FALSE}}
}
        

\seealso{ 
    \code{\link{fitdist}}, \code{\link{bootdist}}.
}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} 
}

%\note{  }

\examples{
# (1) Fit of a lognormal distribution with maximum likelihood estimation
# on toxicity values to estimate (that is called a species sensitivity distribution, 
# SSD, in ecotoxicology, followed by estimation of the 5 percent quantile value of 
# the fitted distribution (that is called the HC5 in ecotoxicology)
#
data(endosulfan)
ATV <-subset(endosulfan,group == "NonArthroInvert")$ATV
fATV <- fitdist(ATV,"lnorm")
qATV5 <- quantiledist(fATV,probs = 0.05,bootstrap=FALSE)
qATV5boot <- quantiledist(fATV,probs = 0.05,bootstrap=TRUE)
qATV5boot

# (2) Estimation of quantiles of 
# the fitted distribution for various probabilities (HCx with various x in ecotoxicology)
#
qATVboot <- quantiledist(fATV,probs = c(0.05,0.1,0.2,0.5),bootstrap=TRUE)
qATVboot

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
