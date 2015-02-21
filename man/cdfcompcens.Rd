\name{cdfcompcens}
\alias{cdfcompcens}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compares various fitted distributions to censored data on a cdf plot}

\description{
 Plots the empirical cumulative distribution (censored continuous data) with theoretical ones
 corresponding to various fitted distributions.
}

\usage{
cdfcompcens(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
    datacol, fitlty, fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", 
    ylegend = NULL, lines01 = FALSE, Turnbull.confint = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{ft}{One \code{"fitdistcens"} object or a list of objects of class \code{"fitdistcens"}.}
\item{xlim}{The \eqn{x}-limits of the plot.}
\item{ylim}{The \eqn{y}-limits of the plot.}
\item{xlogscale}{If \code{TRUE}, uses a logarithmic scale for the \eqn{x}-axis.}
\item{ylogscale}{If \code{TRUE}, uses a logarithmic scale for the \eqn{y}-axis.}
\item{main}{A main title for the plot, see also \code{\link{title}}.}
\item{xlab}{A label for the \eqn{x}-axis, defaults to a description of \code{x}.}
\item{ylab}{A label for the \eqn{y}-axis, defaults to a description of \code{y}.}
\item{datacol}{A specification of the color to be used in plotting data points.}
\item{fitcol}{A (vector of) color(s) to plot fitted distributions. 
    If there are fewer colors than fits they are recycled in the standard fashion.}
\item{fitlty}{A (vector of) line type(s) to plot fitted distributions. 
    If there are fewer colors than fits they are recycled in the standard fashion.
    See also \code{\link{par}}.}
\item{addlegend}{If \code{TRUE}, a legend is added to the plot.}
\item{legendtext}{A character or expression vector of length \eqn{\geq 1} to appear
    in the legend, see also \code{\link{legend}}.}
\item{xlegend, ylegend}{The \eqn{x} and \eqn{y} co-ordinates to be used to position the legend.
          They can be specified by keyword or in any way which is
          accepted by 'xy.coords': see \code{\link{legend}} for details.}
\item{lines01}{A logical to plot two horizontal lines at \code{h=0} and \code{h=1}.}
\item{Turnbull.confint}{ if TRUE, confidence intervals will be added to the Turnbull plot
    of the empirical distribution.}
\item{\dots}{Further graphical arguments passed to graphical functions used in \code{\link{cdfcompcens}}.}
}

\details{
    Empirical and theoretical distributions are plotted in cdf. 
    The EM approach of Turnbull (Turnbull, 1974) is used to compute the overall 
    empirical cdf curve, with confidence intervals if
    \code{Turnbull.confint} is \code{TRUE}, 
    by calls to functions \code{\link{survfit}} and \code{\link{plot.survfit}} from the 
    \code{survival} package.
    By default a legend is added to the plot. Many graphical arguments are optional, 
    dedicated to personalize the plot, and fixed to default values if omitted.
}

\seealso{ 
    \code{\link{plotdistcens}}, \code{\link{survfit.formula}}, \code{\link{legend}} 
    and \code{\link{par}}.
}

\references{ 
Turnbull BW (1974), \emph{Nonparametric estimation of a survivorship function with doubly
censored data}.
Journal of American Statistical Association, 69, 169-173.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr}. 
}

%\note{  }

\examples{
# (1) Plot various distributions fitted to bacterial contamination data
#
data(smokedfish)
Clog10 <- log10(smokedfish)

fitsfn <- fitdistcens(Clog10,"norm")
summary(fitsfn)

fitsfl <- fitdistcens(Clog10,"logis")
summary(fitsfl)

dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q,a,b) exp(-exp((a-q)/b))
qgumbel <- function(p,a,b) a-b*log(-log(p))
fitsfg<-fitdistcens(Clog10,"gumbel",start=list(a=-3,b=3))
summary(fitsfg)

cdfcompcens(list(fitsfn,fitsfl,fitsfg))
cdfcompcens(list(fitsfn,fitsfl,fitsfg),datacol="orange",
    legendtext=c("normal","logistic","Gumbel"),
    main="bacterial contamination fits",
    xlab="bacterial concentration (CFU/g)",ylab="F",
    xlegend = "center",lines01 = TRUE)

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
