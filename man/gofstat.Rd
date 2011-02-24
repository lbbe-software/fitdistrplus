\name{gofstat}
\alias{gofstat}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Goodness-of-fit statistics}

\description{
 Computes goodness-of-fit statistics for a fit of a parametric distribution on 
 non-censored data.  
}

\usage{
gofstat(f, chisqbreaks, meancount, print.test = FALSE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{f}{ An object of class 'fitdist' result of the function \code{fitdist}.}
\item{chisqbreaks}{ A numeric vector defining the breaks of the cells used to compute the chi-squared
    statistic. If omitted, these breaks are automatically computed from the data 
    in order to reach roughly the same number of observations per cell, roughly equal to the argument
    \code{meancount}, or sligthly more if there are some ties. }
\item{meancount}{ The mean number of observations per cell expected for the definition of the breaks 
    of the cells used to compute the chi-squared statistic. This argument will not be taken into 
    account if the breaks are directly defined in the argument \code{chisqbreaks}. If \code{chisqbreaks} 
    and \code{meancount} are both omitted, \code{meancount} is fixed in order to obtain roughly 
    \eqn{(4n)^{2/5}} cells with \eqn{n} the length of the dataset. }
\item{print.test}{ If \code{FALSE}, the results of the tests are computed but not printed  }
}

\details{
     
    Goodness-of-fit statistics are computed. The Chi-squared statistic is computed using cells defined 
    by the argument 
    \code{chisqbreaks} or cells automatically defined from the data in order 
    to reach roughly the same number of observations per cell, roughly equal to the argument
    \code{meancount}, or sligthly more if there are some ties. 
    The choice to define cells from the empirical distribution (data) and not from the
    theoretical distribution was done to enable the comparison of Chi-squared values obtained
    with different distributions fitted on a same dataset.
    If \code{chisqbreaks} and \code{meancount} 
    are both omitted, \code{meancount} is fixed in order to obtain roughly \eqn{(4n)^{2/5}} cells, 
    with \eqn{n} the length of the dataset (Vose, 2000).
    The Chi-squared statistic is not computed if the program fails
    to define enough cells due to a too small dataset. When the Chi-squared statistic is computed, 
    and if the degree of freedom (nb of cells - nb of parameters - 1)  of the corresponding distribution 
    is strictly positive, the p-value of the Chi-squared test is returned.
    
    For the distributions assumed continuous (all but \code{"binom"}, 
    \code{"nbinom"}, \code{"geom"}, \code{"hyper"} and \code{"pois"} for R base distributions), 
    Kolmogorov-Smirnov, Cramer-von Mises and Anderson-Darling  and statistics are also computed, 
    as defined by Stephens (1986). 
    
    An approximate Kolmogorov-Smirnov test is 
    performed by assuming the distribution parameters known. The critical value defined by Stephens (1986) 
    for a completely specified distribution is used to reject or not the 
    distribution at the significance level 0.05. Because of this approximation, the result of the test
    (decision of rejection of the distribution or not) is returned only for datasets with more 
    than 30 observations. Note that this approximate test may be too conservative. 
    
    For datasets with more than 5 observations and for distributions for 
    which the test is described by Stephens (1986) (\code{"norm"}, \code{"lnorm"},
    \code{"exp"}, \code{"cauchy"}, \code{"gamma"}, \code{"logis"} and \code{"weibull"}),
    the Cramer-von Mises and Anderson-darling tests are performed as described by Stephens (1986). 
    Those tests take into 
    account the fact that the parameters are not known but estimated from the data. The result is the 
    decision to reject or not the distribution at the significance level 0.05. Those tests are available
    only for maximum likelihood estimations.
    
    Only recommended statistics are automatically printed, i.e.
    Cramer-von Mises, Anderson-Darling and Kolmogorov statistics for continuous distributions and
    Chi-squared statistics for discrete ones ( \code{"binom"}, 
    \code{"nbinom"}, \code{"geom"}, \code{"hyper"} and \code{"pois"} ). Results of the tests are printed
    only if \code{print.test=TRUE}. Even not printed, all the available results may be found in the list 
    returned by the function. 
    
}

\value{ 
    \code{gof} returns a list with following components,
    \item{ chisq }{ the Chi-squared statistic or \code{NULL} if not computed }
    \item{ chisqbreaks }{ breaks used to define cells in the Chi-squared statistic }
    \item{ chisqpvalue }{ p-value of the Chi-squared statistic or \code{NULL} if not computed}
    \item{ chisqdf }{ degree of freedom of the Chi-squared distribution or \code{NULL} if not computed }
    \item{ chisqtable }{ a table with observed and theoretical counts used for the Chi-squared calculations }
    \item{ cvm }{ the Cramer-von Mises statistic or \code{NULL} if not computed }
    \item{ cvmtest }{ the decision of the Cramer-von Mises test or \code{NULL} if not computed }
    \item{ ad }{ the Anderson-Darling statistic or \code{NULL} if not computed }
    \item{ adtest }{ the decision of the Anderson-Darling test or \code{NULL} if not computed }
    \item{ ks }{ the Kolmogorov-Smirnov statistic or \code{NULL} if not computed }
    \item{ kstest }{ the decision of the Kolmogorov-Smirnov test or \code{NULL} if not computed }
}

\seealso{ 
    \code{\link{fitdist}}.
}

\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 81-155.

Stephens MA (1986) Tests based on edf statistics. In Goodness-of-fit techniques (D'Agostino RB and
Stephens MA, eds), Marcel dekker, New York, pp. 97-194.

Venables WN and Ripley BD (2002) Modern applied statistics with S. Springer, New York, pp. 435-446.

Vose D (2000) Risk analysis, a quantitative guide. John Wiley & Sons Ltd, Chischester, England, pp. 99-143.
}


\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vetagro-sup.fr} and
Christophe Dutang
}

%\note{  }

\examples{

# (1) for a fit of a normal distribution 
#

x1 <- c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
print(f1 <- fitdist(x1,"norm"))
gofstat(f1)
gofstat(f1,print.test=TRUE)

# (2) fit a discrete distribution (Poisson)
#

x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
print(f2<-fitdist(x2,"pois"))
g2 <- gofstat(f2,chisqbreaks=c(0,1),print.test=TRUE)
g2$chisqtable


# (3) comparison of fits of various distributions
#

x3<-rweibull(n=100,shape=2,scale=1)
gofstat(f3a<-fitdist(x3,"weibull"))
gofstat(f3b<-fitdist(x3,"gamma"))
gofstat(f3c<-fitdist(x3,"exp"))

# (4) Use of Chi-squared results in addition to
#     recommended statistics for continuous distributions
#

x4<-rweibull(n=100,shape=2,scale=1)
f4<-fitdist(x4,"weibull")
g4 <-gofstat(f4,meancount=10)
print(g4)

# (5) estimation of the standard deviation of a normal distribution 
# by maximum likelihood with the mean fixed at 10 using the argument fix.arg
#
f1b <- fitdist(x1,"norm",start=list(sd=5),fix.arg=list(mean=10))
gofstat(f1b)


}
\keyword{ distribution }
