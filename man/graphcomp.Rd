\name{graphcomp}
\alias{graphcomp}
\alias{cdfcomp}
\alias{denscomp}
\alias{qqcomp}
\alias{ppcomp}
\title{Graphical comparison of multiple fitted distributions (for non-censored data)}
\description{   
 \code{cdfcomp} plots the empirical cumulative distribution against fitted distribution functions,
 \code{denscomp} plots the histogram against fitted density functions, 
 \code{qqcomp} plots theoretical quantiles against empirical ones,
 \code{ppcomp} plots theoretical probabilities against empirical ones.
 Only \code{cdfcomp} is able to plot fits of a discrete distribution.
}

\usage{
cdfcomp(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
    datapch, datacol, fitlty, fitcol, addlegend = TRUE, legendtext, 
    xlegend = "bottomright", ylegend = NULL, horizontals = TRUE, 
    verticals = FALSE, do.points = TRUE, use.ppoints = TRUE, a.ppoints = 0.5, 
    lines01 = FALSE, discrete, add = FALSE, plotstyle = "graphics", \dots)
    
denscomp(ft, xlim, ylim, probability = TRUE, main, xlab, ylab, datacol, fitlty, 
    fitcol, addlegend = TRUE, legendtext, xlegend = "topright", ylegend = NULL, 
    demp = FALSE, dempcol = "grey", plotstyle = "graphics", \dots)

qqcomp(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
    fitpch, fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", 
    ylegend = NULL, use.ppoints = TRUE, a.ppoints = 0.5, line01 = TRUE, 
    line01col = "black", line01lty = 1, ynoise = TRUE, plotstyle = "graphics", \dots)
    
ppcomp(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
    fitpch, fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", 
    ylegend = NULL, use.ppoints = TRUE, a.ppoints = 0.5, line01 = TRUE, 
    line01col = "black", line01lty = 1, ynoise = TRUE, plotstyle = "graphics", \dots)
}
\arguments{
\item{ft}{One \code{"fitdist"} object or a list of objects of class \code{"fitdist"}.}
\item{xlim}{The \eqn{x}-limits of the plot.}
\item{ylim}{The \eqn{y}-limits of the plot.}
\item{xlogscale}{If \code{TRUE}, uses a logarithmic scale for the \eqn{x}-axis.}
\item{ylogscale}{If \code{TRUE}, uses a logarithmic scale for the \eqn{y}-axis.}
\item{main}{A main title for the plot. See also \code{\link{title}}.}
\item{xlab}{A label for the \eqn{x}-axis, defaults to a description of \code{x}.}
\item{ylab}{A label for the \eqn{y}-axis, defaults to a description of \code{y}.}
\item{datapch}{An integer specifying a symbol to be used in plotting data points. 
  See also \code{\link{par}}.}
\item{datacol}{A specification of the color to be used in plotting data points.
  See also \code{\link{par}}.}
\item{fitcol}{A (vector of) color(s) to plot fitted distributions. 
  If there are fewer colors than fits they are recycled in the standard fashion.
  See also \code{\link{par}}.}
\item{fitlty}{A (vector of) line type(s) to plot fitted distributions/densities. 
  If there are fewer colors than fits they are recycled in the standard fashion.
  See also \code{\link{par}}.}
\item{fitpch}{A (vector of) line type(s) to plot fitted quantiles/probabilities. 
  If there are fewer colors than fits they are recycled in the standard fashion.
  See also \code{\link{par}}.}
\item{addlegend}{If \code{TRUE}, a legend is added to the plot.}
\item{legendtext}{A character or expression vector of length \eqn{\ge 1} to appear
  in the legend. See also \code{\link{legend}}.}
\item{xlegend, ylegend}{The \eqn{x} and \eqn{y} coordinates to be used to position the legend. 
  They can be specified by keyword. 
  If \code{plotstyle = "graphics"}, see \code{\link{xy.coords}} and \code{\link{legend}}.
  If \code{plotstyle = "ggplot"}, the \code{xlegend} keyword must be one of \code{top}, \code{bottom}, 
  \code{left}, or \code{right}. See also \code{guide_legend} in \code{ggplot2}}
\item{horizontals}{If \code{TRUE}, draws horizontal lines for the step empirical 
  cumulative distribution function (ecdf). See also \code{\link{plot.stepfun}}.}
\item{do.points}{If \code{TRUE} (by default), draws points at the x-locations. 
  For large dataset (n > 1e4), \code{do.points} is ignored and no point is drawn.}
\item{verticals}{If \code{TRUE}, draws vertical lines for the empirical cumulative distribution 
  function (ecdf). Only taken into account if \code{horizontals=TRUE}. }
\item{use.ppoints}{If \code{TRUE}, probability points of the empirical distribution
  are defined using function \code{\link{ppoints}} as 
  \code{(1:n - a.ppoints)/(n - 2a.ppoints + 1)}.
  If \code{FALSE}, probability points are simply defined as \code{1:n / n}. This argument is ignored 
  for discrete data.}
\item{a.ppoints}{If \code{use.ppoints=TRUE}, this is passed to the \code{\link{ppoints}} function.}
\item{lines01}{A logical to plot two horizontal lines at \code{h=0} and \code{h=1} for \code{cdfcomp}.}
\item{line01}{A logical to plot an horizontal line \eqn{y=x} for \code{qqcomp} and \code{ppcomp}.}
\item{line01col, line01lty}{Color and line type for \code{line01}. See also \code{\link{par}}.}
\item{demp}{A logical to add the empirical density on the plot, using the 
  \code{\link{density}} function.}
\item{dempcol}{A color for the empirical density in case it is added on the plot (\code{demp=TRUE}).}
\item{ynoise}{A logical to add a small Gaussian noise when plotting empirical 
  quantiles/probabilities for \code{qqcomp} and \code{ppcomp}.}
\item{probability}{A logical to use the probability scale for \code{denscomp}. See also \code{\link{hist}}.}
\item{discrete}{If \code{TRUE}, the distributions are considered discrete.
  If missing, \code{discrete} is passed from the first object of class \code{"fitdist"} of 
  the list \code{ft}.}
\item{add}{If \code{TRUE}, adds to an already existing plot. If \code{FALSE}, starts a new plot.
  This parameter is not available when \code{plotstyle = "ggplot"}.}
\item{plotstyle}{\code{"graphics"} or \code{"ggplot"}. If \code{"graphics"}, the display is built with \code{\link{graphics}} functions. 
  If \code{"ggplot"}, a graphic object output is created with \code{ggplot2} functions.}
  
\item{\dots}{Further graphical arguments passed to graphical functions used in \code{cdfcomp}, \code{denscomp}, \code{ppcomp} 
  and \code{qqcomp} when \code{plotstyle = "graphics"}. 
  When \code{plotstyle = "ggplot"}, these arguments are only used by the histogram plot (\code{\link{hist}}) in the \code{denscomp} function.
  When \code{plotstyle = "ggplot"}, the graphical output can be customized with relevant \code{ggplot2} functions after 
  you store your output.}
}

\details{
  \code{cdfcomp} provides a plot of the empirical distribution and each fitted 
  distribution in cdf, by default using the Hazen's rule
  for the empirical distribution, with probability points defined as 
  \code{(1:n - 0.5)/n}. If \code{discrete} is \code{TRUE}, probability points 
  are always defined as \code{(1:n)/n}. For large dataset (\code{n > 1e4}), no 
  point is drawn but the line for \code{ecdf} is drawn instead.
  Note that when \code{horizontals, verticals and do.points} are \code{FALSE},
  no empirical point is drawn, only the fitted cdf is shown.
    
  \code{denscomp} provides a density plot of each fitted distribution
  with the histogram of the data (data are assumed continuous). 
    
  \code{ppcomp} provides a plot of the probabilities of each fitted distribution
  (\eqn{x}-axis) against the empirical probabilities (\eqn{y}-axis) by default defined as 
  \code{(1:n - 0.5)/n} (data are assumed continuous). 
  For large dataset (\code{n > 1e4}), lines are drawn instead of pointss and customized with the \code{fitpch} parameter.
    
  \code{qqcomp} provides a plot of the quantiles of each theoretical distribution (\eqn{x}-axis)
  against the empirical quantiles of the data (\eqn{y}-axis), by default defining
  probability points as \code{(1:n - 0.5)/n} for theoretical quantile calculation 
  (data are assumed continuous). 
  For large dataset (\code{n > 1e4}), lines are drawn instead of points and customized with the \code{fitpch} parameter.
     
  By default a legend is added to these plots. Many graphical arguments are optional, 
  dedicated to personalize the plots, and fixed to default values if omitted.
}

\seealso{ 
  See \code{\link{CIcdfplot}} and \code{\link{plot}}, \code{\link{legend}}, \code{\link{ppoints}}, 
  \code{\link{plot.stepfun}} for classic plotting functions and \code{\link{plotdist}}.  
}


\references{
Delignette-Muller ML and Dutang C (2015), 
\emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}

\author{ 
Christophe Dutang, Marie-Laure Delignette-Muller and Aurélie Siberchicot.
}

\examples{
# (1) Plot various distributions fitted to serving size data
#
data(groundbeef)
serving <- groundbeef$serving
fitW <- fitdist(serving, "weibull")
fitln <- fitdist(serving, "lnorm")
fitg <- fitdist(serving, "gamma")

cdfcomp(list(fitW, fitln, fitg), horizontals = FALSE)
cdfcomp(list(fitW, fitln, fitg), horizontals = TRUE)
cdfcomp(list(fitW, fitln, fitg), horizontals = TRUE, verticals = TRUE, datacol = "purple")
cdfcomp(list(fitW, fitln, fitg), legendtext = c("Weibull", "lognormal", "gamma"),
  main = "ground beef fits", xlab = "serving sizes (g)",
  ylab = "F", xlim = c(0, 250), xlegend = "center", lines01 = TRUE)
denscomp(list(fitW, fitln, fitg), legendtext = c("Weibull", "lognormal", "gamma"),
  main = "ground beef fits", xlab = "serving sizes (g)", xlim = c(0, 250), xlegend = "topright")
ppcomp(list(fitW, fitln, fitg), legendtext = c("Weibull", "lognormal", "gamma"),
  main = "ground beef fits", xlegend = "bottomright", line01 = TRUE)
qqcomp(list(fitW, fitln, fitg), legendtext = c("Weibull", "lognormal", "gamma"),
  main = "ground beef fits", xlegend = "bottomright", line01 = TRUE,
  xlim = c(0, 300), ylim = c(0, 300), fitpch = 16)


# (2) Plot lognormal distributions fitted by 
# maximum goodness-of-fit estimation
# using various distances (data plotted in log scale)
#
data(endosulfan)
ATV <- subset(endosulfan, group == "NonArthroInvert")$ATV
flnMGEKS <- fitdist(ATV, "lnorm", method = "mge", gof = "KS")
flnMGEAD <- fitdist(ATV, "lnorm", method = "mge", gof = "AD")
flnMGEADL <- fitdist(ATV, "lnorm", method = "mge", gof = "ADL")
flnMGEAD2L <- fitdist(ATV, "lnorm", method = "mge", gof = "AD2L")

cdfcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
  xlogscale = TRUE, main = "fits of a lognormal dist. using various GOF dist.",
  legendtext = c("MGE KS", "MGE AD", "MGE ADL", "MGE AD2L"),
  verticals = TRUE, xlim = c(10, 100000))
qqcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
  main = "fits of a lognormal dist. using various GOF dist.",
  legendtext = c("MGE KS", "MGE AD", "MGE ADL", "MGE AD2L"),
  xlogscale = TRUE, ylogscale = TRUE)
ppcomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
  main = "fits of a lognormal dist. using various GOF dist.",
  legendtext = c("MGE KS", "MGE AD", "MGE ADL", "MGE AD2L"))


# (3) Plot normal and logistic distributions fitted by 
# maximum likelihood estimation
# using various plotting positions in cdf plots
#
data(endosulfan)
log10ATV <-log10(subset(endosulfan, group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV, "norm")
fll <- fitdist(log10ATV, "logis")

# default plot using Hazen plotting position: (1:n - 0.5)/n
cdfcomp(list(fln, fll), legendtext = c("normal", "logistic"), xlab = "log10ATV")

# plot using mean plotting position (named also Gumbel plotting position)
# (1:n)/(n + 1)
cdfcomp(list(fln, fll),legendtext = c("normal", "logistic"), xlab = "log10ATV",
  use.ppoints = TRUE, a.ppoints = 0)
    
# plot using basic plotting position: (1:n)/n
cdfcomp(list(fln, fll),legendtext = c("normal", "logistic"), xlab = "log10ATV",
  use.ppoints = FALSE)


# (4) Comparison of fits of two distributions fitted to discrete data
#
data(toxocara)
number <- toxocara$number
fitp <- fitdist(number, "pois")
fitnb <- fitdist(number, "nbinom")
cdfcomp(list(fitp, fitnb), legendtext = c("Poisson", "negative binomial"))


# (5) Customizing of graphical output and use of ggplot2
#
data(groundbeef)
serving <- groundbeef$serving
fitW <- fitdist(serving, "weibull")
fitln <- fitdist(serving, "lnorm")
fitg <- fitdist(serving, "gamma")
if (requireNamespace ("ggplot2", quietly = TRUE)) {
denscomp(list(fitW, fitln, fitg), plotstyle = "ggplot")
cdfcomp(list(fitW, fitln, fitg), plotstyle = "ggplot")
qqcomp(list(fitW, fitln, fitg), plotstyle = "ggplot")
ppcomp(list(fitW, fitln, fitg), plotstyle = "ggplot")
}

# customizing graphical output with graphics
denscomp(list(fitW, fitln, fitg), legendtext = c("Weibull", "lognormal", "gamma"),
  main = "ground beef fits", xlab = "serving sizes (g)", xlim = c(0, 250),
  xlegend = "topright", addlegend = FALSE)

# customizing graphical output with ggplot2
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  dcomp <- denscomp(list(fitW, fitln, fitg), legendtext = c("Weibull", "lognormal", "gamma"),
    xlab = "serving sizes (g)", xlim = c(0, 250),
    xlegend = "topright", plotstyle = "ggplot", breaks = 20, addlegend = FALSE)
  dcomp + ggplot2::theme_minimal() + ggplot2::ggtitle("Ground beef fits")
}

}

\keyword{ distribution }
