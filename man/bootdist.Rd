\name{bootdist}
\alias{bootdist}
\alias{plot.bootdist}
\alias{print.bootdist}
\alias{summary.bootdist}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Bootstrap simulation of uncertainty for non-censored data}

\description{
 Uses parametric or nonparametric bootstrap resampling in order to simulate uncertainty 
 in the parameters of the distribution fitted to non-censored data.  
}

\usage{
bootdist(f, bootmethod="param", niter=1001)
\method{print}{bootdist}(x,...)
\method{plot}{bootdist}(x,...)
\method{summary}{bootdist}(object,...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{f}{ An object of class 'fitdist' result of the function \code{fitdist}.}
\item{bootmethod}{ A character string coding for the type of resampling : \code{"param"} for a parametric
    resampling and \code{"nonparam"} for a nonparametric resampling of data.}
\item{niter}{ The number of samples drawn by bootstrap.  }
\item{x}{ an object of class 'bootdist'.  }
\item{object}{ an object of class 'bootdist'.  }
\item{\dots}{ further arguments to be passed to generic methods }
}

\details{
    Samples are drawn by parametric bootstrap (resampling from the distribution fitted by \code{fitdist}) or
    non parametric bootstrap (resampling with replacement from the data set). On each bootstrap sample the function 
    \code{mledist}  or the function \code{mmedist} (according to the component \code{f$method} of the 
    object of class 'fitdist') is used to estimate bootstrapped values of parameters. When \code{mledist} fails
    to converge, \code{NA} values are returned. Medians and 2.5 and 97.5 percentiles are computed by removing
    \code{NA} values. The number of iterations for which \code{mledist} converges is printed in the summary with
    the medians and the 95 percent confidence intervals of parameters (2.5 and 97.5 percentiles).
    
    The plot of an object of class 'bootdist' consists in a scatterplot or a matrix of scatterplots
    of the bootstrapped values of parameters.
    It uses the function \code{\link{stripchart}} when the fitted distribution
    is characterized by only one parameter, and the function \code{\link{plot}} in other cases.
    In these last cases, it provides
    a representation of the joint uncertainty distribution of the fitted parameters.
}

\value{ 
    \code{bootdist} returns an object of class 'bootdist', a list with 4 components,
    \item{ estim }{ a data frame containing the boostrapped values of parameters.}
    \item{ converg }{ a vector containing the codes for convergence obtained when using \code{mledist} on 
        each bootstraped data set or \code{NULL} if \code{mmedist} is used.}
    \item{ method }{ A character string coding for the type of resampling : \code{"param"} for a parametric
        resampling and \code{"nonparam"} for a nonparametric resampling of data. }
    \item{ CI }{ bootstrap medians and 95 percent confidence percentile intervals of parameters. }
}
        

\seealso{ 
    \code{\link{fitdist}}, \code{\link{mledist}}, \code{\link{qmedist}} and \code{\link{mmedist}}.
}

\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 181-241.
}


\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vetagro-sup.fr} 
}

%\note{  }

\examples{
# (1) basic fit of a normal distribution with maximum likelihood estimation
# followed by parametric bootstrap
#
x1<-c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
f1<-fitdist(x1,"norm",method="mle")
b1<-bootdist(f1)
print(b1)
plot(b1)
summary(b1)

# (2) non parametric bootstrap
#
b1np<-bootdist(f1,bootmethod="nonparam")
summary(b1np)

# (3) fit of a gamma distribution followed by parametric bootstrap
#
f1b<-fitdist(x1,"gamma",method="mle")
b1b<-bootdist(f1b)
summary(b1b)

# (4) fit of a gamma distribution with control of the optimization
# method, followed by parametric bootstrap
#
f1c <- fitdist(x1,"gamma",optim.method="L-BFGS-B",lower=c(0,0))
b1c <- bootdist(f1c)
summary(b1c)

# (5) estimation of the standard deviation of a normal distribution 
# by maximum likelihood with the mean fixed at 10 using the argument fix.arg
# followed by parametric bootstrap
#
f1d <-fitdist(x1,"norm",start=list(sd=5),fix.arg=list(mean=10))
b1d <- bootdist(f1d)
summary(b1d)
plot(b1d)

# (6) fit of a discrete distribution followed by parametric bootstrap
#
x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
f2<-fitdist(x2,"pois",method="mme")
b2<-bootdist(f2)
plot(b2,pch=16)
summary(b2)

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
