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
datapch,datacol,fitcol,fitlty,xlab,ylab,
main,xlegend,ylegend,...)
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
\item{main}{An overall title for the plot}
\item{xlegend}{The x co-ordinate to be used to position the legend}
\item{ylegend}{The y co-ordinate to be used to position the legend}
\item{\dots}{ further graphical arguments passed to graphical functions used in cdfcomp }
}

\details{
    Empirical and theoretical distributions are simply plotted in cdf. 
    By default a legend is added to the plot. Arguments
    legendtext,datapch,datacol,fitcol,fitlty,xlab,ylab,main,xlegend,ylegend are optional and
    fixed to default values if omitted.
}

\seealso{ 
    \code{\link{plotdist}}.
}



\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vetagro-sup.fr} 
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
cdfcomp(list(fitW,fitln,fitg))
cdfcomp(list(fitW,fitln,fitg),legendtext=c("Weibull","lognormal","gamma"),
    main="ground beef fits",xlab="serving sizes (g)",
    ylab="F")


# (2) Plot in of lognormal distributions fitted by 
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


}

\keyword{ distribution }% at least one, from doc/KEYWORDS
