\name{plotdistcens}
\alias{plotdistcens}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Plot of empirical and theoretical distributions for censored data}

\description{
 Plots an empirical distribution for censored data with a theoretical one if specified.
}

\usage{
plotdistcens(censdata,distr,para,leftNA=-Inf,rightNA=Inf,Turnbull=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{censdata}{ A dataframe of two columns respectively named \code{left} 
    and \code{right}, describing each observed value as an interval.
    The \code{left} column contains either \code{NA} for left censored observations,
     the left bound of the interval for interval censored observations,
     or the observed value for non-censored observations.
     The \code{right} column contains either \code{NA} for right censored observations,
     the right bound of the interval for interval censored observations,
     or the observed value for non-censored observations.  }
\item{distr}{ A character string \code{"name"} naming a distribution, for which the corresponding 
    density function \code{dname} and the corresponding distribution function \code{pname}
    must be defined, or directly the density function.}
\item{para}{ A named list giving the parameters of the named distribution. This argument may be 
    omitted only if \code{distr} is omitted.  }
\item{leftNA}{ the real value of the left bound of left censored observations : \code{-Inf} 
    or a finite value such as \code{0} for positive data for example.}
\item{rightNA}{ the real value of the right bound of right censored observations : \code{Inf}
     or a finite value such as a realistic maximum value.}
\item{Turnbull}{ if TRUE the Turnbull algorithm is used to estimate the cdf curve of the censored data
    and previous arguments \code{leftNA} and \code{rightNA} are not used (see details)}
\item{...}{ further graphical arguments passed to other methods }
}

\details{
    Empirical and, if specified, theoretical distributions are plotted
    in cdf. 
    If \code{Turnbull} is \code{TRUE}, the EM approach of Turnbull (Turnbull, 1974) 
    is used to compute the overall 
    empirical cdf curve with confidence intervals, 
    by calls to functions \code{survfit} and \code{plot.survfit} from the 
    \code{survival} package. 
    
    Else data are reported directly as segments for interval, left and right censored data, 
    and as points for non-censored data. Before plotting, observations are ordered and a rank r
    is associated to each of them. Left censored observations are ordered
    first, by their right bounds. Interval censored and non censored observations
    are then ordered by their mid-points and, at last, right censored observations are
    ordered by their left bounds. If \code{leftNA} (resp. \code{rightNA}) is finite,
    left censored (resp. right censored) observations are considered as interval censored
    observations and ordered by mid-points with non-censored and interval censored data.
    It is sometimes necessary to fix \code{rightNA} or \code{leftNA} to a realistic 
    extreme value, even if not exactly known, to obtain a reasonable global ranking of 
    observations. After ranking, each of the n observations is plotted as a point (one x-value) 
    or a segment (an interval of possible x-values),
    with an y-value equal to r/n, r being the rank of each observation in the global ordering
    previously described. This second method may be interesting but 
    is certainly less rigorous than the Turnbull method
    that should be prefered.
}

\seealso{ 
    \code{\link{plotdist}}, \code{\link{survfit.formula}}.
}

\references{ 
Turnbull BW (1974) Nonparametric estimation of a survivorship function with doubly
censored data.\emph{Journal of American Statistical Association}, \bold{69}, 169-173.
}


\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vetagro-sup.fr} 
}

%\note{  }

\examples{
# (1) Plot of an empirical censored distribution (censored data) as a CDF
# using the default Turnbull method
#
d1<-data.frame(
left=c(1.73,1.51,0.77,1.96,1.96,-1.4,-1.4,NA,-0.11,0.55,
    0.41,2.56,NA,-0.53,0.63,-1.4,-1.4,-1.4,NA,0.13),
right=c(1.73,1.51,0.77,1.96,1.96,0,-0.7,-1.4,-0.11,0.55,
    0.41,2.56,-1.4,-0.53,0.63,0,-0.7,NA,-1.4,0.13))
plotdistcens(d1)
plotdistcens(d1,col="red")

# (2) Add the CDF of a normal distribution 
#
plotdistcens(d1,"norm",para=list(mean=0.12,sd=1.4))

# (3) Basic plot of the same empirical distribution with intervals and points
# defining a realistic maximum value for right censored values
# in the second plot
#
plotdistcens(d1,Turnbull = FALSE)
plotdistcens(d1,rightNA=3, Turnbull = FALSE)

# (4) Plot of the CDF of the same dataset after logarithmic transformation
#   with a lognormal distribution, successively using the two proposed methods
#
d3<-data.frame(left=10^(d1$left),right=10^(d1$right))
plotdistcens(d3,"lnorm",para=list(meanlog=0.27,sdlog=3.3))
plotdistcens(d3,"lnorm",para=list(meanlog=0.27,sdlog=3.3),Turnbull = FALSE, leftNA=0)

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
