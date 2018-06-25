\name{graphcompcens}
\alias{graphcompcens}
\alias{cdfcompcens}\alias{denscompcens}
\alias{qqcompcens}
\alias{ppcompcens}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Graphical comparison of multiple fitted distributions for censored data}

\description{   
 \code{cdfcompcens} plots the empirical cumulative distribution against fitted distribution functions,
 \code{qqcompcens} plots theoretical quantiles against empirical ones,
 \code{ppcompcens} plots theoretical probabilities against empirical ones.
}

\usage{
cdfcompcens(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
    datacol, fillrect, fitlty, fitcol, addlegend = TRUE, legendtext, 
    xlegend = "bottomright", ylegend = NULL, lines01 = FALSE, 
    Turnbull.confint = FALSE, NPMLE.method = "Wang", add = FALSE, \dots)
qqcompcens(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab,
      fillrect, fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", 
      ylegend = NULL, line01 = TRUE, line01col = "black", line01lty = 1, 
      ynoise = TRUE, \dots)
ppcompcens(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
      fillrect, fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", 
      ylegend = NULL, line01 = TRUE, line01col = "black", line01lty = 1, 
      ynoise = TRUE, \dots)
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
\item{fillrect}{A specification of the color to be used for filling rectanges
  of non uniqueness of the empirical cumulative distribution 
  (only used if \code{NPMLE.method} is equal to \code{"Wang"} in \code{cdfcompcens}).
  Fix it to \code{NA} if you do not want to fill the rectangles.}
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
\item{lines01}{A logical to plot two horizontal lines at \code{h=0} and \code{h=1}
  for \code{cdfcompcens}.}
\item{Turnbull.confint}{ if TRUE confidence intervals will be added to the Turnbull plot.
    In that case NPMLE.method is forced to \code{"Turnbull"}}
\item{NPMLE.method}{Two NPMLE technique are provided, \code{"Wang"}, the default one,
    which is implemented in the package npsurv and \code{"Turnbull"}, an older one,
    which is implemented in the package survival (see details). Only \code{"Wang"}
    enables the derivation of a Q-Q plot and a P-P plot.}
\item{add}{logical; if \code{TRUE} add to an already existing plot; if \code{FALSE} start 
  a new plot.}  
\item{line01}{A logical to plot an horizontal line \eqn{y=x} for 
  \code{qqcompcens} and \code{ppcompcens}.}
\item{line01col, line01lty}{Color and line type for \code{line01}. See also \code{\link{par}}.}
\item{ynoise}{A logical to add a small noise when plotting empirical 
  quantiles/probabilities for \code{qqcompcens} and \code{ppcompcens}.}
\item{\dots}{Further graphical arguments passed to graphical functions 
  used in \code{cdfcompcens}, \code{ppcompcens} 
  and \code{qqcompcens}.}
}

\details{
  See details of \code{\link{plotdistcens}} for a detailed description of provided goddness-of-fit plots.
}

\seealso{ 
    \code{\link{plotdistcens}}, \code{\link{survfit.formula}}, \code{\link{legend}}     and \code{\link{par}}.
}

\references{ 
Turnbull BW (1974), \emph{Nonparametric estimation of a survivorship function with doubly
censored data}.
Journal of American Statistical Association, 69, 169-173.

Wang Y (2008), \emph{Dimension-reduced nonparametric maximum likelihood computation
for interval-censored data}. Computational Statistics & Data Analysis, 52, 2388-2402.

Wang Y and Taylor SM (2013), \emph{Efficient computation of nonparametric survival
functions via a hierarchical mixture formulation}. Statistics and Computing, 23, 713-725.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}


\author{ 
Marie-Laure Delignette-Muller and Christophe Dutang.
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

# CDF plot
cdfcompcens(list(fitsfn,fitsfl,fitsfg))
cdfcompcens(list(fitsfn,fitsfl,fitsfg),datacol="orange",fillrect = NA, 
    legendtext=c("normal","logistic","Gumbel"),
    main="bacterial contamination fits",
    xlab="bacterial concentration (CFU/g)",ylab="F",
    xlegend = "bottom",lines01 = TRUE)
# alternative Turnbull plot for the empirical cumulative distribution
# (default plot of the previous versions of the package)
cdfcompcens(list(fitsfn,fitsfl,fitsfg), NPMLE.method = "Turnbull")

# PP plot
ppcompcens(list(fitsfn,fitsfl,fitsfg))
ppcompcens(list(fitsfn,fitsfl,fitsfg), ynoise = FALSE)
par(mfrow = c(2,2))
ppcompcens(fitsfn)
ppcompcens(fitsfl)
ppcompcens(fitsfg)
par(mfrow = c(1,1))

# QQ plot
qqcompcens(list(fitsfn,fitsfl,fitsfg))
qqcompcens(list(fitsfn,fitsfl,fitsfg), ynoise = FALSE)
par(mfrow = c(2,2))
qqcompcens(fitsfn)
qqcompcens(fitsfl)
qqcompcens(fitsfg)
par(mfrow = c(1,1))

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
