\name{smokedfish}
\alias{smokedfish}
\docType{data}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Contamination data of Listeria monocytogenes in smoked fish }

\description{
Contamination data of \emph{Listeria monocytogenes} in smoked fish on the Belgian market in the period
2005 to 2007. 
}

\usage{
data(smokedfish)
}
%- maybe also 'usage' for other objects documented here.

\format{ 
\code{smokedfish} is a data frame with 2 columns named left and right, describing
each observed value of \emph{Listeria monocytogenes} concentration (in CFU/g) as an interval. 
The left column contains either NA for
left censored observations, the left bound of the interval for interval censored
observations, or the observed value for non-censored observations. The right
column contains either NA for right censored observations, the right bound of
the interval for interval censored observations, or the observed value for non-censored
observations.
}

\source{
Busschaert, P., Geereard, A.H., Uyttendaele, M., Van Impe, J.F., 2010. 
Estimating distributions out of qualitative and (semi) quantitative microbiological 
contamination data for use in risk assessment. \emph{International Journal of Food Microbiology}. 
\bold{138}, 260-269.
 }

%\references{  }

\examples{
# (1) load of data
#
data(smokedfish)

# (2) plot of data in CFU/g
#
plotdistcens(smokedfish)

# (3) plot of transformed data in log10[CFU/g]
#
Clog10 <- log10(smokedfish)
plotdistcens(Clog10)

# (4) Fit of a normal distribution to data in log10[CFU/g]
#
fitlog10 <- fitdistcens(Clog10,"norm")
summary(fitlog10)
plot(fitlog10)

}

\keyword{ datasets }% at least one, from doc/KEYWORDS
