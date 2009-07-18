\name{plotdist}
\alias{plotdist}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Plot of empirical and theoretical distributions for non-censored data}

\description{
 Plots an empirical distribution (non-censored data) with a theoretical one if specified.
}

\usage{
plotdist(data,distr,para,breaks="default",discrete=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{data}{ A numeric vector.  }
\item{distr}{ A character string \code{"name"} naming a distribution for which the corresponding 
    density function \code{dname}, the corresponding distribution function \code{pname} and the corresponding
    quantile function \code{qname} must be defined, or directly the density function. 
    This argument may be omitted only if \code{para} is omitted. }
\item{para}{ A named list giving the parameters of the named distribution. This argument may be 
    omitted only if \code{distname} is omitted.  }
\item{breaks}{ If \code{"default"} the histogram is plotted with the function \code{hist} 
    with its default breaks definition. Else \code{breaks} is passed to the function \code{hist}.
    This argument is not taken into account if \code{discrete} is \code{TRUE}. }
\item{discrete}{ If TRUE, the distribution is considered as discrete. This argument is not taken into 
    account if \code{distname} is defined. In this last case, the distribution is automatically assumed
    discrete if and only if \code{distname} is \code{"binom"}, \code{"nbinom"},\code{"geom"},
    \code{"hyper"} or \code{"pois"}.  }
\item{...}{ further graphical arguments passed to other methods }
}

\details{
    Empirical and, if specified, theoretical distributions are plotted
    in density and in cdf. For continuous distributions, the function \code{\link{hist}} is used with its default
    breaks definition if \code{breaks} is \code{"default"} or passing \code{breaks} as an argument if it differs
    from \code{"default"}. For continuous distribution and when a theoretical distribution is specified
    by both arguments \code{distname} and \code{para}, Q-Q plot and P-P plot are also plotted (Cullen and Fray, 1999).
    The function \code{\link{ppoints}} is used for the Q-Q plot, to generate the set of probabilities at 
    which to evaluate the inverse distribution.
}

\seealso{ 
    \code{\link{descdist}}, \code{\link{hist}}, \code{\link{plotdistcens}} and \code{\link{ppoints}}.
}

\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 81-159.
}

\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vet-lyon.fr} 
}

%\note{  }

\examples{
x1<-c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
plotdist(x1)

x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
plotdist(x2,discrete=TRUE)
plotdist(x2,"pois",para=list(lambda=mean(x2)))

xn<-rnorm(n=100,mean=10,sd=5)
plotdist(xn,"norm",para=list(mean=mean(xn),sd=sd(xn)))
}

\keyword{ distribution }% at least one, from doc/KEYWORDS