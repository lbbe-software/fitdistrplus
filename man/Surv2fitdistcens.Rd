\name{Surv2fitdistcens}
\alias{Surv2fitdistcens}
\title{Data handling based on Surv() to prepare for fitdistcens()}

\description{
 Provide a function to prepare a data.frame needed by fitdistcens
}

\usage{
Surv2fitdistcens(time, time2, event,
                             type=c('right', 'left', 'interval'))
    
}

\arguments{
\item{time}{for right censored data, this is the follow up time. 
  For interval data, the first argument is the starting time for the interval.}
\item{event}{The status indicator, normally \code{0}=alive, \code{1}=dead. Other choices 
  are \code{TRUE/FALSE} (\code{TRUE} = death) or \code{1/2} (\code{2}=death). 
  For interval censored data, the status indicator is \code{0}=right censored, 
  \code{1}=event at time, \code{2}=left censored, \code{3}=interval censored. 
  For factor data, assume that it has only two levels with the second level
  coding death.}
\item{time2}{ending time of the interval for interval censored or counting 
  process data only. Intervals are assumed to be open on the left and closed 
  on the right, (start, end]. For counting process data, event indicates 
  whether an event occurred at the end of the interval.}
\item{type}{character string specifying the type of censoring. Possible 
  values are \code{"right"}, \code{"left"}, \code{"interval"}.}

}

\details{
    \code{Surv2fitdistcens} makes a \code{data.frame} with two columns
    respectively named \code{left} and \code{right}, describing each observed 
    value as an interval. The \code{left} column contains either \code{NA} 
    for left-censored observations, the left bound of the interval for 
    interval-censored observations, or the observed value for non-censored observations.
    The right column contains either \code{NA} for right-censored observations, 
    the right bound of the interval for interval censored observations, 
    or the observed value for non-censored observations..

}

\value{ 
    \code{fitdistcens} returns a data.frame with two columns 
    respectively named \code{left} and \code{right}.
}


\seealso{ 
    See \code{\link{fitdistrplus}} for an overview of the package.
    See \code{\link{fitdistcens}} for fitting of univariate distributions to censored data
    and \code{\link{fremale}} for the full dataset used in examples below.
    See \code{\link[survival:Surv]{Surv}} for survival objects which use the same arguments.
    
 
}

\references{ 

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}


\author{ 
Christophe Dutang and Marie-Laure Delignette-Muller.
}

%\note{  }

\examples{
# (1) randomized fictive data for male survival - right-censored
#

origdata <- rbind(
c(   43.01164, 55.00890,     0),
c(   36.37782, 47.17864,     0),
c(   33.10335, 34.51882,     0),
c(   71.00890, 81.15537,     1),
c(   80.89254, 81.91376,     1),
c(   67.81656, 78.48323,     1),
c(   73.98494, 76.92539,     1),
c(   71.21971, 79.21971,     1),
c(   68.96920, 79.92608,     1),
c(   53.19370, 54.80082,     1))
colnames(origdata) <- c("AgeIn", "AgeOut", "Death")

Surv2fitdistcens(origdata[,"AgeIn"], origdata[,"AgeOut"], origdata[,"Death"])

# (2) randomized fictive data for male survival - interval-censored
#
origdata[6:7, "Death"] <- 2
origdata[8:10, "Death"] <- 3

Surv2fitdistcens(origdata[,"AgeIn"], origdata[,"AgeOut"], origdata[,"Death"],
  type = "interval")
}

\keyword{ distribution }% at least one, from doc/KEYWORDS
