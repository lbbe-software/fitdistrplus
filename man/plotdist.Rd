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
    omitted only if \code{distr} is omitted.  }
\item{breaks}{ If \code{"default"} the histogram is plotted with the function \code{hist} 
    with its default breaks definition. Else \code{breaks} is passed to the function \code{hist}.
    This argument is not taken into account if \code{discrete} is \code{TRUE}. }
\item{discrete}{ If TRUE, the distribution is considered as discrete. This argument is not taken into 
    account if \code{distr} is defined. In this last case, the distribution is automatically assumed
    discrete if and only if \code{distr} is \code{"binom"}, \code{"nbinom"},\code{"geom"},
    \code{"hyper"} or \code{"pois"}.  }
\item{\dots}{ further graphical arguments passed to graphical functions used in plotdist }
}

\details{
    Empirical and, if specified, theoretical distributions are plotted
    in density and in cdf. For continuous distributions, the function \code{\link{hist}} is used with its default
    breaks definition if \code{breaks} is \code{"default"} or passing \code{breaks} as an argument if it differs
    from \code{"default"}. For continuous distribution and when a theoretical distribution is specified
    by both arguments \code{distname} and \code{para}, Q-Q plot
    (plot of the quantiles of the theoretical fitted distribution (x-axis) against the empirical quantiles of the data)
    and P-P plot (i.e. for each value of the data set, plot of the cumulative density function of the fitted distribution
    (x-axis) against the empirical cumulative density function (y-axis)) are also given (Cullen and Frey, 1999).
    The function \code{\link{ppoints}} (with default parameter for argument a)
    is used for the Q-Q plot, to generate the set of probabilities at 
    which to evaluate the inverse distribution.
    NOTE THAT FROM VERSION 0.4-3, \code{\link{ppoints}} is also used for P-P plot and cdf plot for continuous data.
    To personalize the four plots proposed for continuous data, for example to change the plotting position, we recommend 
    the use of functions \code{\link{cdfcomp}}, \code{denscomp}, \code{QQcomp} and \code{PPcomp} 
    (the three last ones will be soon available).
}

\seealso{ 
    \code{\link{cdfcomp}}, \code{\link{descdist}}, \code{\link{hist}},  \code{\link{plot}}, \code{\link{plotdistcens}} and \code{\link{ppoints}}.
}


\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 81-159.
}

\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vetagro-sup.fr} 
}

%\note{  }

\examples{
# (1) Plot of an empirical distribution with changing 
# of default line types for CDF and colors
#
x1 <- c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
plotdist(x1)
plotdist(x1,col="red",type="b",pch=16)
plotdist(x1,type="s")

# (2) Plot of a discrete distribution against data
#
x2 <- c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
plotdist(x2,discrete=TRUE)
plotdist(x2,"pois",para=list(lambda=mean(x2)))
plotdist(x2,"pois",para=list(lambda=mean(x2)),col="red",lwd="2")

# (3) Plot of a continuous distribution against data
#
xn <- rnorm(n=100,mean=10,sd=5)
plotdist(xn,"norm",para=list(mean=mean(xn),sd=sd(xn)))
plotdist(xn,"norm",para=list(mean=mean(xn),sd=sd(xn)),pch=16,col="green")

# (4) Plot of serving size data
#
data(groundbeef)
plotdist(groundbeef$serving)
}

\keyword{ distribution }% at least one, from doc/KEYWORDS
