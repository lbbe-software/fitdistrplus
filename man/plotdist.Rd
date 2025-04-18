\name{plotdist}
\alias{plotdist}

\title{ Plot of empirical and theoretical distributions for non-censored data}

\description{
 Plots an empirical distribution (non-censored data) with a theoretical one if specified.
}

\usage{
plotdist(data, distr, para, histo = TRUE, breaks = "default", 
  demp = FALSE, discrete, \dots)
}

\arguments{
\item{data}{ A numeric vector.  }
\item{distr}{ A character string \code{"name"} naming a distribution for which the corresponding 
    density function \code{dname}, the corresponding distribution function \code{pname} and the corresponding
    quantile function \code{qname} must be defined, or directly the density function. 
    This argument may be omitted only if \code{para} is omitted. }
\item{para}{ A named list giving the parameters of the named distribution. This argument may be 
    omitted only if \code{distr} is omitted.  }
\item{histo}{A logical to plot the histogram using the 
   \code{\link{hist}} function.}
\item{breaks}{ If \code{"default"} the histogram is plotted with the function \code{hist} 
    with its default breaks definition. Else \code{breaks} is passed to the function \code{hist}.
    This argument is not taken into account if \code{discrete} is \code{TRUE}. }
\item{demp}{A logical to plot the empirical density on the first plot 
  (alone or superimposed on the histogram depending of the value of the argument \code{histo}) 
  using the \code{\link{density}} function.}
\item{discrete}{ If TRUE, the distribution is considered as discrete.
	If both \code{distr} and \code{discrete} are missing, \code{discrete} is set
	to \code{FALSE}. If \code{discrete} is missing but not  \code{distr}, 
	\code{discrete} is set to \code{TRUE} when \code{distr} belongs to 
	\code{"binom"}, \code{"nbinom"},\code{"geom"},
    \code{"hyper"} or \code{"pois"}.  }
\item{\dots}{ further graphical arguments passed to graphical functions used in plotdist.}
}

\details{
    Empirical and, if specified, theoretical distributions are plotted
    in density and in cdf. For the plot in density, the user can use the arguments
    \code{histo} and \code{demp} to specify if he wants the histogram using the function 
    \code{\link{hist}}, the density plot using the function \code{\link{density}}, or both
    (at least one of the two arguments must be put to \code{"TRUE"}).
    For continuous distributions, the function \code{\link{hist}} is used with its default
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
    the use of functions \code{\link{cdfcomp}}, \code{\link{denscomp}}, \code{\link{qqcomp}} and \code{\link{ppcomp}}.
}

\seealso{ 
    See \code{\link{graphcomp}}, \code{\link{descdist}}, \code{\link{hist}},  
    \code{\link{plot}}, \code{\link{plotdistcens}} and \code{\link{ppoints}}.
    
    Please visit the \href{https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html}{Frequently Asked Questions}.
}


\references{ 
Cullen AC and Frey HC (1999), \emph{Probabilistic techniques in exposure assessment}.
Plenum Press, USA, pp. 81-155.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34, \doi{https://doi.org/10.18637/jss.v064.i04}.

}

\author{ 
Marie-Laure Delignette-Muller and Christophe Dutang.
}

\examples{
# (1) Plot of an empirical distribution with changing 
# of default line types for CDF and colors
# and optionally adding a density line
#
set.seed(1234)
x1 <- rnorm(n=30)
plotdist(x1)
plotdist(x1,demp = TRUE)
plotdist(x1,histo = FALSE, demp = TRUE)
plotdist(x1, col="blue", type="b", pch=16)
plotdist(x1, type="s")

# (2) Plot of a discrete distribution against data
#
set.seed(1234)
x2 <- rpois(n=30, lambda = 2)
plotdist(x2, discrete=TRUE)
plotdist(x2, "pois", para=list(lambda = mean(x2)))
plotdist(x2, "pois", para=list(lambda = mean(x2)), lwd="2")

# (3) Plot of a continuous distribution against data
#
xn <- rnorm(n=100, mean=10, sd=5)
plotdist(xn, "norm", para=list(mean=mean(xn), sd=sd(xn)))
plotdist(xn, "norm", para=list(mean=mean(xn), sd=sd(xn)), pch=16)
plotdist(xn, "norm", para=list(mean=mean(xn), sd=sd(xn)), demp = TRUE)
plotdist(xn, "norm", para=list(mean=mean(xn), sd=sd(xn)), 
histo = FALSE, demp = TRUE)

# (4) Plot of serving size data
#
data(groundbeef)
plotdist(groundbeef$serving, type="s")

# (5) Plot of numbers of parasites with a Poisson distribution 
data(toxocara)
number <- toxocara$number
plotdist(number, discrete = TRUE)
plotdist(number,"pois",para=list(lambda=mean(number)))

}

\keyword{ distribution }
