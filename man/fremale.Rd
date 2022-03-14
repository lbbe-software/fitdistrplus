\name{fremale}
\alias{fremale}
\docType{data}

\title{Fictive survival dataset of a french Male population}

\description{100 male individuals randomly taken from \code{frefictivetable} in \code{CASdatasets} package}

\usage{data(fremale)}

\format{ 
\code{fremale} is a data frame with 3 columns names \code{AgeIn}, \code{AgeOut}
respectively for entry age and exit age; \code{Death} a binary dummy: 1 indicating
the death of the individual; 0 a censored observation.
}

\references{See full dataset \code{frefictivetable} of \code{CASdatasets} at \url{http://dutangc.perso.math.cnrs.fr/RRepository/}}

\examples{
# (1) load of data
#
data(fremale)
summary(fremale)
}

\keyword{datasets}
