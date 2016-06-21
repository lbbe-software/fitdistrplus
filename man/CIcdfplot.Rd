\name{CIcdfplot}
\alias{CIcdfplot}
\title{Empirical cumulative distribution function with pointwise confidence intervals on probabilities or on quantiles}
\description{   
 \code{cdfband} plots the empirical cumulative distribution function with the bootstraped pointwise confidence intervals on probabilities of on quantiles.
}

\usage{
CIcdfplot(b, CI.output, CI.type = "two.sided", CI.level = 0.95, CI.col = "pink", CI.lty = 2, 
  CI.fill = NULL, CI.only = FALSE, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, 
  xlab, ylab, datapch, datacol, fitlty, fitcol, addlegend = FALSE, legendtext, 
  xlegend = "bottomright", ylegend = NULL, horizontals = TRUE, verticals = FALSE, 
  do.points = TRUE, use.ppoints = TRUE, a.ppoints = 0.5, lines01 = FALSE, \dots)
}
\arguments{
\item{b}{One \code{"bootdist"} object.}
\item{CI.output}{ The quantity on which (bootstraped) bootstraped confidence intervals are computed: 
    either \code{"probability"} or \code{"quantile"}).}
\item{CI.type}{ Type of confidence intervals : either \code{"two.sided"} or one-sided 
    intervals (\code{"less"} or \code{"greater"}).}
\item{CI.level}{ The confidence level.}
\item{CI.col}{the color of the confidence intervals.}    
\item{CI.lty}{the line type of the confidence intervals.}    
\item{CI.fill}{a color to fill the confidence area. Default is \code{NULL} corresponding to no filling.}    
\item{CI.only}{A logical whether to plot empirical and fitted distribution functions 
    or only the confidence intervals. Default to \code{FALSE}.}

\item{xlim}{The \eqn{x}-limits of the plot.}
\item{ylim}{The \eqn{y}-limits of the plot.}
\item{xlogscale}{If \code{TRUE}, uses a logarithmic scale for the \eqn{x}-axis.}
\item{ylogscale}{If \code{TRUE}, uses a logarithmic scale for the \eqn{y}-axis.}
\item{main}{A main title for the plot, see also \code{\link{title}}.}
\item{xlab}{A label for the \eqn{x}-axis, defaults to a description of \code{x}.}
\item{ylab}{A label for the \eqn{y}-axis, defaults to a description of \code{y}.}

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
    \code{CIcdfplot} provides a plot of the empirical distribution with bootstraped pointwise confidence        intervals on probabilities (y values) or on quantiles (x values). 
    Each interval is computed by evaluating the quantity of interest (probability
    associated to an x value or quantile associated to an y value) using all the
    bootstraped values of parameters to get a bootstraped sample
    of the quantity of interest and then by calculating percentiles on this sample to get a
    confidence interval (classically 2.5 and 97.5 percentiles for a 95 percent 
    confidence level).
    If \code{CI.fill != NULL}, then the whole confidence area is filled by the color \code{CI.fill}
    thanks to the function \code{polygon}, otherwise only borders are drawn thanks to the function
    \code{matline}. Further graphical arguments can be passed to these functions using 
    the three dots arguments \code{\dots}.
}

\seealso{ 
    See also \code{\link{cdfcomp}}, \code{\link{bootdist}} and \code{\link{quantile}}.  
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
# We choose a low number of bootstrap replicates in order to satisfy CRAN running times
# constraint.
# For practical applications, we recommend to use at least niter=501 or niter=1001.


# (1) Fit of an exponential distribution
#

set.seed(123)
s1 <- rexp(50, 1)
f1 <- fitdist(s1, "exp")
b1 <- bootdist(f1, niter= 51) 

# plot 95 percent bilateral confidence intervals on y values (probabilities)
CIcdfplot(b1, CI.level= 95/100, CI.output = "probability")

# plot of the previous intervals as a band 
CIcdfplot(b1, CI.level= 95/100, CI.output = "probability", 
  CI.fill = "pink", CI.col = "red")

# plot of the previous intervals as a band without empirical and fitted dist. functions
CIcdfplot(b1, CI.level= 95/100, CI.output = "probability", CI.only = TRUE,
  CI.fill = "pink", CI.col = "red")

# plot 95 percent bilateral confidence intervals on x values (quantiles)
# with coloring the band
CIcdfplot(b1, CI.level= 95/100, CI.output = "quantile",
  CI.fill = "pink", CI.col = "red")


# plot 95 percent unilateral confidence intervals on quantiles
CIcdfplot(b1, CI.level = 95/100, CI.output = "quant", CI.type = "less", 
  CI.fill = "grey80", CI.col = "black")
CIcdfplot(b1, CI.level= 95/100, CI.output = "quant", CI.type = "greater",
  CI.fill = "grey80", CI.col = "black")

# (2) Fit of a normal distribution on acute toxicity log-transformed values of 
# endosulfan for nonarthropod invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5, 10 and 20 percent quantile  
# values of the fitted distribution, which are called the 5, 10, 20 percent hazardous 
# concentrations (HC5, HC10, HC20) in ecotoxicology, with their
# confidence intervals, from a small number of bootstrap 
# iterations to satisfy CRAN running times constraint and plot of the band
# representing pointwise confidence intervals on any quantiles (any HCx values)
# For practical applications, we recommend to use at least niter=501 or niter=1001.
#

data(endosulfan)
ATV <- subset(endosulfan, group == "NonArthroInvert")$ATV
log10ATV <- log10(subset(endosulfan, group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV, "norm")
bln <- bootdist(fln, bootmethod ="param", niter=101)
quantile(bln, probs = c(0.05, 0.1, 0.2))
CIcdfplot(bln, CI.output = "quantile", CI.fill = "lightblue", CI.col = "blue")

}

\keyword{ distribution }
