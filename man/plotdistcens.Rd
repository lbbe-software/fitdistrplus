\name{plotdistcens}
\alias{plotdistcens}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Plot of empirical and theoretical distributions for censored data}

\description{
 Plots an empirical distribution for censored data with a theoretical one if specified.
}

\usage{
plotdistcens(censdata, distr, para, leftNA = -Inf, rightNA = Inf,NPMLE = TRUE,
    Turnbull = NULL, Turnbull.confint = FALSE, NPMLE.method = "Wang", \dots)
}

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
\item{NPMLE}{ if TRUE an NPMLE (nonparametric maximum likelihood estimate) technique is 
    used to estimate the cdf curve of the censored data
    and previous arguments \code{leftNA} and \code{rightNA} are not used (see details)}
\item{Turnbull}{ Turnbull is a deprecated argument. You should now use \code{NPMLE} to choose
    to use an NPMLE technique (or not) and \code{NPMLE.method} to choose the method.}
\item{Turnbull.confint}{ if TRUE confidence intervals will be added to the Turnbull plot.
    In that case NPMLE.method is forced to \code{"Turnbull"}}
\item{NPMLE.method}{Two NPMLE technique are provided, \code{"Wang"}, the default one,
    which is implemented in the package npsurv and \code{"Turnbull"}, an older one,
    which is implemented in the package survival (see details). Only \code{"Wang"}
    enables the derivation of a Q-Q plot and a P-P plot.}
\item{\dots}{ further graphical arguments passed to other methods }
}

\details{
     If \code{NPMLE} is \code{TRUE}, and \code{NPMLE.method} is \code{"Wang"},
    empirical and, if specified, theoretical distributions are plotted
    in cdf using either the constrained Newton method (Wang, 2008)
    or the hierarchical constrained Newton method (Wang, 2013) 
    to compute the overall 
    empirical cdf curve, 
    by calls to functions \code{\link{npsurv}} from the 
    \code{npsurv} package. In this plot, grey rectangles represent areas
    where the empirical distribution function is not unique. In cases
    where a theoretical distribution is specified, two goodness-of-fit plots
    are also provided, a Q-Q plot (plot of the quantiles of the theoretical fitted 
    distribution (x-axis) against the empirical quantiles of the data) and a P-P plot
    (i.e. for each value of the data set, plot of the cumulative density function 
    of the fitted distribution (x-axis) against the empirical cumulative density function 
    (y-axis)). Grey rectangles in a Q-Q plot or a P-P plot also represent areas of 
    non uniqueness of empirical quantiles or probabilities, directly derived from 
    non uniqueness areas of the empirical cumulative distribution.

    If \code{NPMLE} is \code{TRUE}, and \code{NPMLE.method} is \code{"Turnbull"},
    empirical and, if specified, theoretical distributions are plotted
    in cdf using the EM approach of Turnbull (Turnbull, 1974) 
    to compute the overall 
    empirical cdf curve, with confidence intervals if \code{Turnbull.confint} is \code{TRUE}, 
    by calls to functions \code{\link{survfit}} and \code{\link{plot.survfit}} from the 
    \code{survival} package. 
    
    If \code{NPMLE} is \code{FALSE} empirical and, if specified, theoretical distributions 
    are plotted in cdf, with data directly reported as segments for interval, 
    left and right censored data, 
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
Turnbull BW (1974), \emph{Nonparametric estimation of a survivorship function with doubly
censored data}. Journal of American Statistical Association, 69, 169-173.

Wang Y (2008), \emph{Dimension-reduced nonparametric maximum likelihood computation
for interval-censored data}. Computational Statistics & Data Analysis, 52, 2388-2402.

Wang Y and Taylor SM (2013), \emph{Efficient computation of nonparametric survival
functions via a hierarchical mixture formulation}. Statistics and Computing, 23, 713-725.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}


\author{ 
Marie-Laure Delignette-Muller and Christophe Dutang.
}

%\note{  }

\examples{
# (1) Plot of an empirical censored distribution (censored data) as a CDF
# using the default Turnbull method
#
data(smokedfish)
d1 <- as.data.frame(log10(smokedfish))
plotdistcens(d1)

# (2) Add the CDF of a normal distribution 
#
plotdistcens(d1,"norm",para=list(mean=-1.6,sd=1.5))

# (3) Various plots of the same empirical distribution 
#
# default Wang plot
plotdistcens(d1, NPMLE = TRUE, NPMLE.method = "Wang")
# Turnbull plot
plotdistcens(d1, NPMLE = TRUE, NPMLE.method = "Turnbull")
# Turnbull plot with confidence intervals
plotdistcens(d1, NPMLE = TRUE, NPMLE.method = "Turnbull", Turnbull.confint = TRUE)
# with intervals and points
plotdistcens(d1,rightNA=3, NPMLE = FALSE)
# with intervals and points
# defining a minimum value for left censored values
plotdistcens(d1,leftNA=-3, NPMLE = FALSE)

}

\keyword{ distribution }% at least one, from doc/KEYWORDS
