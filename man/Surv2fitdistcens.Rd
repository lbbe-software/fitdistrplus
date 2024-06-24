\name{Surv2fitdistcens}
\alias{Surv2fitdistcens}
\title{Handling of data formated as in the survival package for use in fitdistcens()}

\description{
 Provide a function to prepare a data frame needed by fitdistcens() from 
 data classically coded when using the Surv() function of the survival
 package
}

\usage{
Surv2fitdistcens(time, time2, event,
                      type = c('right', 'left', 'interval', 'interval2'))
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
\item{time2}{ending time of the interval for interval censored. 
  Intervals are assumed to be open on the left and closed 
  on the right, (start, end].}
\item{type}{character string specifying the type of censoring. Possible 
  values are \code{"right"}, \code{"left"}, \code{"interval"}, \code{"interval2"}.}
}

\details{
    \code{Surv2fitdistcens} makes a \code{data.frame} with two columns
    respectively named \code{left} and \code{right}, describing each observed 
    value as an interval as required in fitdistcens(): 
    the \code{left} column contains either \code{NA} 
    for left-censored observations, the left bound of the interval for 
    interval-censored observations, or the observed value for non-censored observations.
    The right column contains either \code{NA} for right-censored observations, 
    the right bound of the interval for interval censored observations, 
    or the observed value for non-censored observations.
}

\value{ 
    \code{Surv2fitdistcens} returns a data.frame with two columns 
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
Journal of Statistical Software, 64(4), 1-34, \doi{https://doi.org/10.18637/jss.v064.i04}.
}


\author{ 
Christophe Dutang and Marie-Laure Delignette-Muller.
}

\examples{
# (1) randomized fictive survival data - right-censored
#
origdata <- data.frame(rbind(
c(   43.01, 55.00,     0),
c(   36.37, 47.17,     0),
c(   33.10, 34.51,     0),
c(   71.00, 81.15,     1),
c(   80.89, 81.91,     1),
c(   67.81, 78.48,     1),
c(   73.98, 76.92,     1),
c(   53.19, 54.80,     1)))
colnames(origdata) <- c("AgeIn", "AgeOut", "Death")

# add of follow-up time (for type = "right" in Surv())
origdata$followuptime <- origdata$AgeOut - origdata$AgeIn
origdata

### use of default survival type "right"
# in Surv()
survival::Surv(time = origdata$followuptime, event = origdata$Death, type = "right")
# for fitdistcens()
Surv2fitdistcens(origdata$followuptime, event = origdata$Death, type = "right")

# use of survival type "interval" 
# in Surv()
survival::Surv(time = origdata$followuptime, time2 = origdata$followuptime, 
          event = origdata$Death, type = "interval")
# for fitdistcens()
Surv2fitdistcens(time = origdata$followuptime, time2 = origdata$followuptime, 
          event = origdata$Death, type = "interval") 

# use of survival type "interval2" 
origdata$survivalt1 <- origdata$followuptime
origdata$survivalt2 <- origdata$survivalt1
origdata$survivalt2[1:3] <- Inf
origdata
survival::Surv(time = origdata$survivalt1, time2 = origdata$survivalt2, 
type = "interval2")
Surv2fitdistcens(origdata$survivalt1, time2 = origdata$survivalt2, 
                type = "interval2")


# (2) Other examples with various left, right and interval censored values
#
# with left censored data
(d1 <- data.frame(time = c(2, 5, 3, 7), ind = c(0, 1, 1, 1)))
survival::Surv(time = d1$time, event = d1$ind, type = "left")
Surv2fitdistcens(time = d1$time, event = d1$ind, type = "left")

(d1bis <- data.frame(t1 = c(2, 5, 3, 7), t2 = c(2, 5, 3, 7), 
  censtype = c(2, 1, 1, 1)))
survival::Surv(time = d1bis$t1, time2 = d1bis$t2, 
  event = d1bis$censtype, type = "interval")
Surv2fitdistcens(time = d1bis$t1, time2 = d1bis$t2, 
  event = d1bis$censtype, type = "interval")

# with interval, left and right censored data
(d2 <- data.frame(t1 = c(-Inf, 2, 3, 4, 3, 7), t2 = c(2, 5, 3, 7, 8, Inf)))
survival::Surv(time = d2$t1, time2 = d2$t2, type = "interval2")
Surv2fitdistcens(time = d2$t1, time2 = d2$t2, type = "interval2")

(d2bis <- data.frame(t1 = c(2, 2, 3, 4, 3, 7), t2 = c(2, 5, 3, 7, 8, 7), 
  censtype = c(2,3,1,3,3,0)))
survival::Surv(time = d2bis$t1, time2 = d2bis$t2, 
  event = d2bis$censtype, type = "interval")
Surv2fitdistcens(time = d2bis$t1, time2 = d2bis$t2, 
  event = d2bis$censtype, type = "interval")
}

\keyword{distribution}