\name{CIcdfplot}
\alias{CIcdfplot}
\title{Empirical cumulative distribution function with confidence interval}
\description{   
 \code{cdfband} plots the empirical cumulative distribution function with the bootstraped confidence interval.
}

\usage{
CIcdfplot(b, CI.output, CI.type = "two.sided", CI.level = 0.95, CI.col="red", CI.lty=2, 
  CI.fill=FALSE, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
  dataplot=TRUE, datapch, datacol, fitlty, fitcol, addlegend = FALSE, legendtext, 
  xlegend = "bottomright", ylegend = NULL, horizontals = TRUE, verticals = FALSE, 
  do.points = TRUE, use.ppoints = TRUE, a.ppoints = 0.5, lines01 = FALSE, \dots)
}
\arguments{
\item{b}{One \code{"bootdist"} object.}
\item{CI.output}{ The quantity on which (bootstraped) quantiles are computed: 
    either \code{"probability"} or \code{"quantile"}).}
\item{CI.type}{ Type of confidence intervals : either \code{"two.sided"} or one-sided 
    intervals (\code{"less"} or \code{"greater"}).}
\item{CI.level}{ The confidence level.}
\item{CI.col}{the color of the confidence interval.}    
\item{CI.lty}{the line type of the confidence interval.}    
\item{CI.fill}{a logical to fill the confidence area. Default is \code{FALSE}.}    

\item{xlim}{The \eqn{x}-limits of the plot.}
\item{ylim}{The \eqn{y}-limits of the plot.}
\item{xlogscale}{If \code{TRUE}, uses a logarithmic scale for the \eqn{x}-axis.}
\item{ylogscale}{If \code{TRUE}, uses a logarithmic scale for the \eqn{y}-axis.}
\item{main}{A main title for the plot, see also \code{\link{title}}.}
\item{xlab}{A label for the \eqn{x}-axis, defaults to a description of \code{x}.}
\item{ylab}{A label for the \eqn{y}-axis, defaults to a description of \code{y}.}

\item{dataplot}{A logical whether to plot \code{ecdf} or not. Default to \code{TRUE}.}
\item{datapch}{An integer specifying a symbol to be used in plotting data points, 
    see also \code{\link{points}}.}
\item{datacol}{A specification of the color to be used in plotting data points.}
\item{fitcol}{A (vector of) color(s) to plot fitted distributions. 
    If there are fewer colors than fits they are recycled in the standard fashion.}
\item{fitlty}{A (vector of) line type(s) to plot fitted distributions/densities. 
    If there are fewer colors than fits they are recycled in the standard fashion.
    See also \code{\link{par}}.}
\item{addlegend}{If \code{TRUE}, a legend is added to the plot.}
\item{legendtext}{A character or expression vector of length \eqn{\geq 1} to appear
    in the legend, see also \code{\link{legend}}.}
\item{xlegend, ylegend}{The \eqn{x} and \eqn{y} co-ordinates to be used to position 
  the legend. They can be specified by keyword or in any way which is
          accepted by 'xy.coords': see \code{\link{legend}} for details.}
\item{horizontals}{If \code{TRUE}, draws horizontal lines for the step empirical 
  cdf function. See also \code{\link{plot.stepfun}}.}
\item{do.points}{logical; if \code{TRUE}, also draw points at the x-locations. 
  Default is true.}
\item{verticals}{If \code{TRUE}, draws also vertical lines for the empirical cdf function. 
    Only taken into account if \code{horizontals=TRUE}. }
\item{use.ppoints}{If \code{TRUE}, probability points of the empirical distribution
    are defined using function \code{\link{ppoints}} as 
    \code{(1:n - a.ppoints)/(n - 2a.ppoints + 1)}.
     If \code{FALSE}, 
    probability points are simply defined as \code{1:n / n}. This argument is ignored 
    for discrete data.}
\item{a.ppoints}{If \code{use.ppoints=TRUE}, this is passed to function \code{\link{ppoints}}.}
\item{lines01}{A logical to plot two horizontal lines at \code{h=0} and \code{h=1} for \code{cdfcomp}.}
\item{\dots}{Further graphical arguments passed to \code{matlines} or \code{polygon},
  respectively when \code{CI.fill=FALSE} and \code{CI.fill=TRUE}.}
}

\details{
    \code{cdfband} provides a plot of the empirical distribution with bootstraped confidence interval.
    If \code{CI.fill=TRUE}, then the whole confidence area is filled by the color \code{CI.col}
    thanks to the function \code{polygon}, otherwise only borders are drawn thanks to the function
    \code{matline}. Further graphical arguments can be passed to these functions using 
    the three dots arguments \code{\dots}.
}

\seealso{ 
    See also \code{\link{cdfcomp}}, \code{\link{bootdist}}.  
}


\references{
Delignette-Muller ML and Dutang C (2015), 
\emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}

\author{ 
Christophe Dutang and Marie-Laure Delignette-Muller.
}

\examples{


# (1) Fit of an exponential distribution
#

set.seed(123)
s1 <- rexp(50, 1)
f1 <- fitdist(s1, "exp")
b1 <- bootdist(f1, niter=11)

# plot 90 percent confidence area
CIcdfplot(b1, CI.level=90/100, CI.output="probability")
CIcdfplot(b1, CI.level=95/100, CI.output="prob", CI.type = "l", 
  CI.fill=TRUE, CI.col="grey80")
CIcdfplot(b1, CI.level=99/100, CI.output="quant", CI.type = "greater")

}

\keyword{ distribution }
