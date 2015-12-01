\name{logLik}
\alias{llsurface}
\alias{llcurve}
\alias{logLik}

\title{3D/2D graphics of likelihood surfaces or 2D plots of likelihood curves}
\description{   
 \code{llsurface} plots the likelihood surface for distributions with two or more parameters,
 \code{llcurve} plots the likelihood curve for distributions with one or more parameters.
}

\usage{
llsurface(plot.min, plot.max, plot.arg, plot.np=50, fix.arg=NULL, obs, distr, 
    loglik=TRUE, plot.type=c("contour", "filled.contour", "persp", "image"), 
    enhance=TRUE, \dots)    

llcurve(plot.min, plot.max, plot.arg, plot.np=50, fix.arg=NULL, obs, distr, 
    loglik=TRUE, plot.type=c("line", "points"), enhance=TRUE, \dots)

}
\arguments{
\item{plot.min}{a two-element vector with lower plotting bounds for 
  \code{llsurface}, only one element for \code{llcurve}.}
\item{plot.max}{a two-element vector with upper plotting bounds for 
  \code{llsurface}, only one element for \code{llcurve}.}
\item{plot.arg}{a two-element vector with the two parameters that will 
  vary for \code{llsurface}, only one element for \code{llcurve}.}
\item{plot.np}{the number of points.}
\item{fix.arg}{a named list with fixed value of other parameters.}
\item{obs}{the vector of observations.}
\item{distr}{the name of distributions.}
\item{loglik}{a logical to plot log-likelihood or likelihood function.}
\item{plot.type}{a character string among \code{"contour"}, 
  \code{"filled.contour"}, \code{"persp"}, \code{"image"} for \code{llsurface},
  a character string among \code{"line"}, \code{"points"} for \code{llcurve}.
}
\item{enhance}{a logical to get an enhance plot.}
\item{\dots}{Further graphical arguments passed to graphical functions.}

}

\details{


\code{llsurface} plots the likelihood surface for distributions with two 
varying parameters and other parameters fixed.
When \code{plot.type="contour"}, \code{\link{contour}} (2D-plot) is used.
When \code{plot.type="filled.contour"}, \code{\link{filled.contour}} 
(2D-plot) is used.
When \code{plot.type="persp"}, \code{\link{persp}} (3D-plot) is used.
When \code{plot.type="image"}, \code{\link{image}} (2D-plot) is used.

\code{llcurve} plots the likelihood curve for distributions with one 
varying parameter and other parameters fixed.
When \code{plot.type="line"}, \code{\link{plot}(type="l")} (2D-plot) is used.
When \code{plot.type="points"}, \code{\link{plot}(type="p")} (2D-plot) is used.

}

\seealso{ 
    See \code{\link{plot}}, \code{\link{contour}}, \code{\link{filled.contour}},
    \code{\link{persp}}, \code{\link{image}} for classic plotting functions.  
}


\references{
Delignette-Muller ML and Dutang C (2015), 
\emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}

\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} 
and Christophe Dutang.
}

\examples{
# (1) likelihood curves
#

n <- 100
set.seed(1234)
x <- rexp(n)

llcurve(0, 4, "rate", distr="exp", obs=x, 
  main="log-likelihood for exponential distribution")
abline(v=1, lty=2)

llcurve(0, 4, "shape", distr="gamma", fix.arg=list(rate=1), obs=x, 
  main="log-likelihood for gamma distribution")
abline(v=1, lty=2)


# (2) likelihood surface
# 

x <- rnorm(n, 0, 1)

par(mfrow=c(1, 3))
llsurface(plot.min=c(-1, 0.5), plot.max=c(1, 3/2), plot.arg=c("mean", "sd"),
                obs=x, distr="norm", plot.type="persp", theta=60, phi=30,
                main="log-likelihood for normal distribution")
llsurface(plot.min=c(-1, 0.5), plot.max=c(1, 3/2), plot.arg=c("mean", "sd"),
                obs=x, distr="norm", plot.type="contour",
                main="log-likelihood for normal distribution")
points(0, 1, pch="+", col="red")
llsurface(plot.min=c(-1, 0.5), plot.max=c(1, 3/2), plot.arg=c("mean", "sd"),
                obs=x, distr="norm", plot.type="image",
                main="log-likelihood for normal distribution")
points(0, 1, pch="+")
                



}

\keyword{ distribution }
