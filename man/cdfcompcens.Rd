\name{cdfcompcens}
\alias{cdfcompcens}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compares various fitted distributions to censored data on a cdf plot}

\description{
 Plots the empirical cumulative distribution (censored continuous data) with theoretical ones
 corresponding to various fitted distributions.
}

\usage{
cdfcompcens(ft,xlogscale=FALSE,addlegend=TRUE,legendtext,datacol,
fitcol,fitlty,xlab,ylab,main,xlegend = "bottomright",
ylegend = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{ft}{ A list of object of class 'fitdist'.  }
\item{xlogscale}{ If TRUE the cdf plot uses a log scale for data (x axis) }
\item{addlegend}{ If TRUE a legend is added to the plot to link colors and types of lines
    to the fitted distributions }
\item{legendtext}{An expression vector of length the number of fits to compare, so 
    the length of list ft }
\item{datacol}{A specification of the color to be 
    used in plotting data points (ecdf)}
\item{fitcol}{A vector of length the number of fits to compare specifying the
     colors of lines used to plot fitted cumulative distributions}
\item{fitlty}{A vector of length the number of fits to compare specifying the
     types of lines used to plot fitted cumulative distributions}
\item{xlab}{A title for the x axis}
\item{ylab}{A title for the y axis}
\item{main}{An overall title for the plot}
\item{xlegend,ylegend}{the x and y co-ordinates to be used to position the legend.
          They can be specified by keyword or in any way which is
          accepted by 'xy.coords': See \code{\link{legend}} for details.}
\item{\dots}{ further graphical arguments passed to graphical functions used in cdfcompcens }
}

\details{
    Empirical and theoretical distributions are simply plotted in cdf. 
    The EM approach of Turnbull (Turnbull, 1974) is used to compute the overall 
    empirical cdf curve with confidence intervals, 
    by calls to functions \code{survfit} and \code{plot.survfit} from the 
    \code{survival} package.
    By default a legend is added to the plot. Arguments
    legendtext,datacol,fitcol,fitlty,xlab,ylab,main,xlegend,ylegend are optional and
    fixed to default values if omitted.
}

\seealso{ 
    \code{\link{plotdistcens}}, \code{\link{survfit.formula}}, \code{\link{legend}}.
}

\references{ 
Turnbull BW (1974) Nonparametric estimation of a survivorship function with doubly
censored data.\emph{Journal of American Statistical Association}, \bold{69}, 169-173.
}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} 
}

%\note{  }

\examples{
# (1) Plot various distributions fitted to bacterial contamination data
#
data(smokedfish)
fitsfn <- fitdistcens(smokedfish,"norm")
summary(fitsfn)

fitsfl <- fitdistcens(smokedfish,"logis")
summary(fitsfl)

dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q,a,b) exp(-exp((a-q)/b))
qgumbel <- function(p,a,b) a-b*log(-log(p))
fitsfg<-fitdistcens(smokedfish,"gumbel",start=list(a=-3,b=3))
summary(fitsfg)

cdfcompcens(list(fitsfn,fitsfl,fitsfg))
cdfcompcens(list(fitsfn,fitsfl,fitsfg),datacol="orange",
    legendtext=c("normal","logistic","Gumbel"),
    main="bacterial contamination fits",
    xlab="bacterial concentration (CFU/g)",ylab="F")
cdfcompcens(list(fitsfn,fitsfl,fitsfg),datacol="orange",
    legendtext=c("normal","logistic","Gumbel"),
    main="bacterial contamination fits",
    xlab="bacterial concentration (CFU/g)",ylab="F",
    xlegend = "center")

# (2) Same plot in x logscale
#
cdfcompcens(list(fitsfn,fitsfl,fitsfg),xlog = TRUE)


}

\keyword{ distribution }% at least one, from doc/KEYWORDS
