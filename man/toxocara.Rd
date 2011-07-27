\name{toxocara}
\alias{toxocara}
\docType{data}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Parasite abundance in insular feral cats }

\description{
Toxocara cati abundance in feral cats living on Kerguelen island. 
}

\usage{
data(toxocara)
}
%- maybe also 'usage' for other objects documented here.

\format{ 
\code{toxocara} is a data frame with 1 column (number: number of parasites in digestive tract)
}

\source{ 
Fromont, E., Morvilliers, L., Artois, M., Pontier, D. 2001. Parasite richness and abundance in insular 
and mainland feral cats. \emph{Parasitology}, \bold{123}, 143-151.
 }

%\references{  }

\examples{
# (1) load of data
#
data(toxocara)

# (2) description and plot of data
#
number <- toxacara$number
descdist(number,discrete=TRUE,boot=1000)
plotdist(number,discrete=TRUE)

# (3) fit of a Poisson distribution to data
#
fitp <- fitdist(number,"pois")
summary(fitp)
plot(fitp)

# (4) fit of a negative binomial distribution to data
#
fitnb <- fitdist(number,"nbinom")
summary(fitnb)
plot(fitnb)

}

\keyword{ datasets }% at least one, from doc/KEYWORDS
