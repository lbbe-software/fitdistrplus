\name{cdfcomp}
\alias{cdfcomp}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compares various fitted distributions to data on a cdf plot}

\description{
 Plots the empirical cumulative distribution (non-censored continuous data) with theoretical ones
 corresponding to various fitted distributions.
}

\usage{
cdfcomp(ft,xlogscale=FALSE,addlegend=TRUE,legendtext,
datapch,datacol,fitcol,fitlty,xlab,ylab,xlim,
main,xlegend = "bottomright",ylegend = NULL,horizontals = TRUE, 
verticals = FALSE, use.ppoints = TRUE, a.ppoints = 0.5,...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{ft}{ A list of object of class 'fitdist'.  }
\item{xlogscale}{ If TRUE the cdf plot uses a log scale for data (x axis) }
\item{addlegend}{ If TRUE a legend is added to the plot to link colors and types of lines
    to the fitted distributions }
\item{legendtext}{An expression vector of length the number of fits to compare, so 
    the length of list ft }
\item{datapch}{An integer specifying a symbol to be 
    used in plotting data points (ecdf)}
\item{datacol}{A specification of the color to be 
    used in plotting data points (ecdf)}
\item{fitcol}{A vector of length the number of fits to compare specifying the
     colors of lines used to plot fitted cumulative distributions}
\item{fitlty}{A vector of length the number of fits to compare specifying the
     types of lines used to plot fitted cumulative distributions}
\item{xlab}{A title for the x axis}
\item{ylab}{A title for the y axis}
\item{xlim}{numeric vectors of length 2, giving the x coordinates range}
\item{main}{An overall title for the plot}
\item{xlegend,ylegend}{the x and y co-ordinates to be used to position the legend.
          They can be specified by keyword or in any way which is
          accepted by 'xy.coords': See \code{\link{legend}} for details.}
\item{horizontals}{If TRUE, draws horizontal lines for step empirical cdf function}
\item{verticals}{If TRUE, draws also vertical lines for step empirical cdf function. Only taken
into account if \code{horizontals} is TRUE }
\item{use.ppoints}{If TRUE, probability points are defined as (1:n - a.ppoints)/(n - 2a.ppoints + 1) 
    for ordinates of the cdf plot of data, using function \code{\link{ppoints}}. If FALSE, 
    probability points are simply defined as 1:n / n.}
\item{a.ppoints}{If \code{use.ppoints} is TRUE, this argument is passed to function \code{\link{ppoints}}.}
\item{\dots}{ further graphical arguments passed to graphical functions used in cdfcomp }
}

\details{
    Empirical and theoretical distributions are simply plotted in cdf by default using the Hazen's rule
    for the empirical distribution, 
    with probability points defined as (1:n - 0.5)/n. 
    By default a legend is added to the plot. Arguments
    legendtext,datapch,datacol,fitcol,fitlty,xlab,ylab,xlim,main,xlegend,ylegend are optional and
    fixed to default values if omitted.
}

\seealso{ 
    \code{\link{plotdist}}, \code{\link{legend}}, \code{\link{ppoints}}.
}



\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} 
}

%\note{  }

\examples{
# (1) Plot various distributions fitted to serving size data
#
data(groundbeef)
serving <- groundbeef$serving
fitW <- fitdist(serving,"weibull")
fitln <- fitdist(serving,"lnorm")
fitg <- fitdist(serving,"gamma")
cdfcomp(list(fitW,fitln,fitg),horizontals = FALSE)
cdfcomp(list(fitW,fitln,fitg),horizontals = TRUE)
cdfcomp(list(fitW,fitln,fitg),horizontals = TRUE, verticals=TRUE,datacol="orange")
cdfcomp(list(fitW,fitln,fitg),legendtext=c("Weibull","lognormal","gamma"),
    main="ground beef fits",xlab="serving sizes (g)",
    ylab="F",xlim = c(0,250))
cdfcomp(list(fitW,fitln,fitg),legendtext=c("Weibull","lognormal","gamma"),
    main="ground beef fits",xlab="serving sizes (g)",
    ylab="F",xlim = c(0,250),xlegend = "topleft")

# (2) Plot lognormal distributions fitted by 
# maximum goodness-of-fit estimation
# using various distances (data plotted in log scale)
#
data(endosulfan)
ATV <-subset(endosulfan,group == "NonArthroInvert")$ATV
flnMGEKS <- fitdist(ATV,"lnorm",method="mge",gof="KS")
flnMGEAD <- fitdist(ATV,"lnorm",method="mge",gof="AD")
flnMGEADL <- fitdist(ATV,"lnorm",method="mge",gof="ADL")
flnMGEAD2L <- fitdist(ATV,"lnorm",method="mge",gof="AD2L")
cdfcomp(list(flnMGEKS,flnMGEAD,flnMGEADL,flnMGEAD2L),
xlogscale=TRUE,main="fits of a lognormal dist. using various GOF dist.",
legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"))
cdfcomp(list(flnMGEKS,flnMGEAD,flnMGEADL,flnMGEAD2L),
xlogscale=TRUE,verticals=TRUE,xlim=c(10,100000))

# (3) Plot normal and logistic distributions fitted by 
# maximum likelihood estimation
# using various plotting positions
#
data(endosulfan)
log10ATV <-log10(subset(endosulfan,group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV,"norm")
fll <- fitdist(log10ATV,"logis")
# default plot using Hazen plotting position: (1:n - 0.5)/n
cdfcomp(list(fln,fll),legendtext=c("normal","logistic"),xlab="log10ATV")
# plot using mean plotting position (named also Gumbel plotting position)
# (1:n)/(n + 1)
cdfcomp(list(fln,fll),legendtext=c("normal","logistic"),xlab="log10ATV",
    use.ppoints = TRUE, a.ppoints = 0)
# plot using basic plotting position: (1:n)/n
cdfcomp(list(fln,fll),legendtext=c("normal","logistic"),xlab="log10ATV",
    use.ppoints = FALSE)

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
