\name{bootdistcens}
\alias{bootdistcens}
\alias{plot.bootdistcens}
\alias{print.bootdistcens}
\alias{summary.bootdistcens}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Bootstrap simulation of uncertainty for censored data}

\description{
 Uses nonparametric bootstrap resampling in order to simulate uncertainty 
 in the parameters of the distribution fitted to censored data.  
}

\usage{
bootdistcens(f, niter=1001)
\method{print}{bootdistcens}(x,...)
\method{plot}{bootdistcens}(x,...)
\method{summary}{bootdistcens}(object,...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{f}{ An object of class 'fitdistcens' result of the function \code{fitdistcens}.}
\item{niter}{ The number of samples drawn by bootstrap.  }
\item{x}{ an object of class 'bootdistcens'.  }
\item{object}{ an object of class 'bootdistcens'.  }
\item{...}{ further arguments to be passed to generic methods }
}

\details{
    Samples are drawn by 
    non parametric bootstrap (resampling with replacement from the data set). On each bootstrap sample the function 
    \code{mledist} is used to estimate bootstrapped values of parameters. When \code{mledist} fails
    to converge, \code{NA} values are returned. Medians and 2.5 and 97.5 percentiles are computed by removing
    \code{NA} values. The number of iterations for which \code{mledist} converges is printed in the summary with
    the medians and the 95 percent confidence intervals of parameters (2.5 and 97.5 percentiles).
    
    The plot of an object of class 'bootdistcens' consists in a scatterplot or a matrix of scatterplots
    of the bootstrapped values of parameters.
    It uses the function \code{\link{stripchart}} when the fitted distribution
    is characterized by only one parameter, and the function \code{\link{plot}} in other cases.
    In these last cases, it provides
    a representation of the joint uncertainty distribution of the fitted parameters.
}

\value{ 
    \code{bootdistcens} returns an object of class 'bootdistcens', a list with 3 components,
    \item{ estim }{ a data frame containing the boostrapped values of parameters.}
    \item{ converg }{ a vector containing the codes for convergence obtained when using \code{mledist} on 
        each bootstraped data set or \code{NULL} if \code{mmedist} is used.}
    \item{ CI }{ bootstrap medians and 95 percent confidence percentile intervals of parameters. }
}
        

\seealso{ 
    \code{\link{fitdistcens}} and \code{\link{mledist}}.
}

\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 181-241.
}


\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vet-lyon.fr} 
}

%\note{  }

\examples{
# (1) fit of a normal distribution followed by nonparametric bootstrap
#
d1<-data.frame(
left=c(1.73,1.51,0.77,1.96,1.96,-1.4,-1.4,NA,-0.11,0.55,
    0.41,2.56,NA,-0.53,0.63,-1.4,-1.4,-1.4,NA,0.13),
right=c(1.73,1.51,0.77,1.96,1.96,0,-0.7,-1.4,-0.11,0.55,
    0.41,2.56,-1.4,-0.53,0.63,0,-0.7,NA,-1.4,0.13))
f1<-fitdistcens(d1, "norm")
b1<-bootdistcens(f1)
b1
summary(b1)
plot(b1)

# (2) fit of a gamma distribution followed by nonparametric bootstrap
#
d3<-data.frame(left=10^(d1$left),right=10^(d1$right))
f3 <- fitdistcens(d3,"gamma")
b3 <- bootdistcens(f3)
summary(b3)
plot(b3)

# (3) fit of a gamma distribution followed by nonparametric bootstrap
# with control of the optimization method
#
f3BFGS <- fitdistcens(d3,"gamma",optim.method="L-BFGS-B",lower=c(0,0))
b3BFGS <- bootdistcens(f3BFGS)
summary(b3BFGS)
plot(b3BFGS)

# (4) estimation of the standard deviation of a normal distribution 
# by maximum likelihood with the mean fixed at 0.1 using the argument fix.arg
# followed by nonparametric bootstrap
#
f1b <- fitdistcens(d1, "norm", start=list(sd=1.5),fix.arg=list(mean=0.1))
b1b<-bootdistcens(f1b)
summary(b1b)
plot(b1b)

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
