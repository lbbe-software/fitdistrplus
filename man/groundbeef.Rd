\name{groundbeef}
\alias{groundbeef}
\docType{data}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Ground beef serving size data set }

\description{
Serving sizes collected in a French survey, for ground beef patties consumed by children under 5 years old. 
}

\usage{
data(groundbeef)
}
%- maybe also 'usage' for other objects documented here.

\format{ 
\code{groundbeef} is a data frame with 1 column (serving: serving sizes in grams)
}

\source{ 
Delignette-Muller, M.L., Cornu, M. 2008. Quantitative risk assessment for Escherichia coli O157:H7 in 
frozen ground beef patties consumed by young children in French households. \emph{International Journal 
of Food Microbiology}, \bold{128}, 158-164.
 }

%\references{  }

\examples{
# (1) load of data
#
data(groundbeef)

# (2) description and plot of data
#
serving <- groundbeef$serving
descdist(serving)
plotdist(serving)

# (3) fit of a Weibull distribution to data
#
fitW <- fitdist(serving,"weibull")
summary(fitW)
plot(fitW)
gofstat(fitW)

}

\keyword{ datasets }% at least one, from doc/KEYWORDS
