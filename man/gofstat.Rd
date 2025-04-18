\name{gofstat}
\alias{gofstat}
\alias{print.gofstat.fitdist}
\alias{print.gofstat.fitdistcens}
\title{ Goodness-of-fit statistics}

\description{
 Computes goodness-of-fit statistics for parametric distributions fitted to a same  
 censored or non-censored data set.  
}

\usage{
gofstat(f, chisqbreaks, meancount, discrete, fitnames=NULL) 
	
\method{print}{gofstat.fitdist}(x, \dots)
\method{print}{gofstat.fitdistcens}(x, \dots)
	
}

\arguments{
\item{f}{ An object of class \code{"fitdist"} (or \code{"fitdistcens"} ), 
output of the function \code{fitdist()} (resp. \code{"fitdist()"}), or
	a list of \code{"fitdist"} objects, or a list of \code{"fitdistcens"} objects.}
\item{chisqbreaks}{Only usable for non censored data, a numeric vector defining the breaks of the cells used to compute the chi-squared
    statistic. If omitted, these breaks are automatically computed from the data 
    in order to reach roughly the same number of observations per cell, roughly equal to the argument
    \code{meancount}, or sligthly more if there are some ties. }
\item{meancount}{Only usable for non censored data, the mean number of observations per cell expected for the definition of the breaks 
    of the cells used to compute the chi-squared statistic. This argument will not be taken into 
    account if the breaks are directly defined in the argument \code{chisqbreaks}. If \code{chisqbreaks} 
    and \code{meancount} are both omitted, \code{meancount} is fixed in order to obtain roughly 
    \eqn{(4n)^{2/5}} cells with \eqn{n} the length of the dataset. }
\item{discrete}{If \code{TRUE}, only the Chi-squared statistic and information criteria are computed.
	If missing, \code{discrete} is passed from the first object of class \code{"fitdist"} of the list \code{f}.
	For censored data this argument is ignored, as censored data are considered continuous.}
\item{fitnames}{A vector defining the names of the fits.}
\item{x}{An object of class \code{"gofstat.fitdist"} or \code{"gofstat.fitdistcens"}.  }
\item{\dots}{Further arguments to be passed to generic functions.}
}

\details{
    For any type of data (censored or not), information criteria are calculated.
    For non censored data, added Goodness-of-fit statistics are computed as described below. 
    
    The Chi-squared statistic is computed using cells defined 
    by the argument 
    \code{chisqbreaks} or cells automatically defined from data, in order 
    to reach roughly the same number of observations per cell, roughly equal to the argument
    \code{meancount}, or sligthly more if there are some ties. 
    The choice to define cells from the empirical distribution (data), and not from the
    theoretical distribution, was done to enable the comparison of Chi-squared values obtained
    with different distributions fitted on a same data set.
    If \code{chisqbreaks} and \code{meancount} 
    are both omitted, \code{meancount} is fixed in order to obtain roughly \eqn{(4n)^{2/5}} cells, 
    with \eqn{n} the length of the data set (Vose, 2000).
    The Chi-squared statistic is not computed if the program fails
    to define enough cells due to a too small dataset. When the Chi-squared statistic is computed, 
    and if the degree of freedom (nb of cells - nb of parameters - 1)  of the corresponding distribution 
    is strictly positive, the p-value of the Chi-squared test is returned.
    
    For continuous distributions, Kolmogorov-Smirnov, Cramer-von Mises and 
	Anderson-Darling  and statistics are also computed, as defined by Stephens (1986). 
    
    An approximate Kolmogorov-Smirnov test is 
    performed by assuming the distribution parameters known. The critical value defined by Stephens (1986) 
    for a completely specified distribution is used to reject or not the 
    distribution at the significance level 0.05. Because of this approximation, the result of the test
    (decision of rejection of the distribution or not) is returned only for data sets with more 
    than 30 observations. Note that this approximate test may be too conservative. 
    
    For data sets with more than 5 observations and for distributions for 
    which the test is described by Stephens (1986) for maximum likelihood estimations 
    (\code{"exp"}, \code{"cauchy"}, \code{"gamma"} and \code{"weibull"}),
    the Cramer-von Mises and Anderson-darling tests are performed as described by Stephens (1986). 
    Those tests take into 
    account the fact that the parameters are not known but estimated from the data by maximum likelihood.
    The result is the 
    decision to reject or not the distribution at the significance level 0.05. Those tests are available
    only for maximum likelihood estimations.
    
    Only recommended statistics are automatically printed, i.e.
    Cramer-von Mises, Anderson-Darling and Kolmogorov statistics for continuous distributions and
    Chi-squared statistics for discrete ones ( \code{"binom"}, 
    \code{"nbinom"}, \code{"geom"}, \code{"hyper"} and \code{"pois"} ). 
	
	Results of the tests are not printed but stored in the output of the function. 
	    
}

\value{ 
    \code{gofstat()} returns an object of class \code{"gofstat.fitdist"} or \code{"gofstat.fitdistcens"} with following components or a sublist of them (only aic, bic and nbfit for censored data) ,
    \item{ chisq }{ a named vector with the Chi-squared statistics or \code{NULL} if not computed}
    \item{ chisqbreaks }{ common breaks used to define cells in the Chi-squared statistic }
    \item{ chisqpvalue }{ a named vector with the p-values of the Chi-squared statistic 
    or \code{NULL} if not computed}
    \item{ chisqdf }{ a named vector with the degrees of freedom of the Chi-squared distribution
    or \code{NULL} if not computed }
    \item{ chisqtable }{ a table with observed and theoretical counts used for the Chi-squared calculations }
    \item{ cvm }{ a named vector of the Cramer-von Mises statistics or \code{"not computed"} if not computed }
    \item{ cvmtest }{ a named vector of the decisions of the Cramer-von Mises test 
    or \code{"not computed"} if not computed }
    \item{ ad }{ a named vector with the Anderson-Darling statistics or 
    \code{"not computed"} if not computed }
    \item{ adtest }{ a named vector with the decisions of the Anderson-Darling test 
    or \code{"not computed"} if not computed }
    \item{ ks }{ a named vector with the Kolmogorov-Smirnov statistic or 
    \code{"not computed"} if not computed }
    \item{ kstest }{ a named vector with the decisions of the Kolmogorov-Smirnov test 
    or \code{"not computed"} if not computed }
	\item{aic}{a named vector with the values of the Akaike's Information Criterion.}
	\item{bic}{a named vector with the values of the Bayesian Information Criterion.}
	\item{discrete}{the input argument or the automatic definition by the function from  the first
    object of class \code{"fitdist"} of the list in input.}
	\item{nbfit}{Number of fits in argument.}
}

\seealso{ 
    See \code{\link{fitdist}}.
    
    Please visit the \href{https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html}{Frequently Asked Questions}.
}

\references{ 
Cullen AC and Frey HC (1999), \emph{Probabilistic techniques in exposure assessment}.
Plenum Press, USA, pp. 81-155.

Stephens MA (1986), \emph{Tests based on edf statistics}. In Goodness-of-fit techniques (D'Agostino RB and
Stephens MA, eds), Marcel Dekker, New York, pp. 97-194.

Venables WN and Ripley BD (2002), \emph{Modern applied statistics with S}.
Springer, New York, pp. 435-446, \doi{10.1007/978-0-387-21706-2}.

Vose D (2000), \emph{Risk analysis, a quantitative guide}.
John Wiley & Sons Ltd, Chischester, England, pp. 99-143.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34, \doi{https://doi.org/10.18637/jss.v064.i04}.

}


\author{ 
Marie-Laure Delignette-Muller and Christophe Dutang.
}


\examples{

# (1) fit of two distributions to the serving size data
# by maximum likelihood estimation
# and comparison of goodness-of-fit statistics
#

data(groundbeef)
serving <- groundbeef$serving
(fitg <- fitdist(serving, "gamma"))
gofstat(fitg)
(fitln <- fitdist(serving, "lnorm"))
gofstat(fitln)

gofstat(list(fitg, fitln))


# (2) fit of two discrete distributions to toxocara data
# and comparison of goodness-of-fit statistics
#

data(toxocara)
number <- toxocara$number

fitp <- fitdist(number,"pois")
summary(fitp)
plot(fitp)

fitnb <- fitdist(number,"nbinom")
summary(fitnb)
plot(fitnb)

gofstat(list(fitp, fitnb),fitnames = c("Poisson","negbin"))

# (3) Get Chi-squared results in addition to
#     recommended statistics for continuous distributions
#

set.seed(1234)
x4 <- rweibull(n=1000,shape=2,scale=1)
# fit of the good distribution
f4 <- fitdist(x4,"weibull")
plot(f4)

# fit of a bad distribution
f4b <- fitdist(x4,"cauchy")
plot(f4b)

(g <- gofstat(list(f4,f4b),fitnames=c("Weibull", "Cauchy")))
g$chisq
g$chisqdf
g$chisqpvalue
g$chisqtable

# and by defining the breaks
(g <- gofstat(list(f4,f4b), 
chisqbreaks = seq(from = min(x4), to = max(x4), length.out = 10), fitnames=c("Weibull", "Cauchy")))
g$chisq
g$chisqdf
g$chisqpvalue
g$chisqtable

# (4) fit of two distributions on acute toxicity values 
# of fluazinam (in decimal logarithm) for
# macroinvertebrates and zooplancton
# and comparison of goodness-of-fit statistics
#

data(fluazinam)
log10EC50 <-log10(fluazinam)
(fln <- fitdistcens(log10EC50,"norm"))
plot(fln)
gofstat(fln)
(fll <- fitdistcens(log10EC50,"logis"))
plot(fll)
gofstat(fll)

gofstat(list(fll, fln), fitnames = c("loglogistic", "lognormal"))


}
\keyword{ distribution }
