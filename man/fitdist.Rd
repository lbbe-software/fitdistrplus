\name{fitdist}
\alias{fitdist}
\alias{plot.fitdist}
\alias{print.fitdist}
\alias{summary.fitdist}
\alias{fitdistrplus}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fitting of univariate distributions to non-censored data and goodness-of-fit statistics}

\description{
 Fits a univariate distribution to non-censored data by maximum likelihood or matching moments,
  and computes goodness-of-fit statistics.  
}

\usage{
fitdist(data, distr, method="mle", start, chisqbreaks, meancount)
\method{print}{fitdist}(x,...)
\method{plot}{fitdist}(x,breaks="default",...)
\method{summary}{fitdist}(object,...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{data}{ A numeric vector.  }
\item{distr}{ A character string \code{"name"} naming a distribution for which the corresponding 
    density function \code{dname}, the corresponding distribution function \code{pname} and the corresponding
    quantile function \code{qname} must be defined, or directly the density function.}
\item{method}{ A character string coding for the fitting method: 
    \code{"mle"} for 'maximum likelihood' and \code{"mom"} for 'matching moments'.}
\item{start}{ A named list giving the initial values of parameters of the named distribution. 
    This argument will not be taken into account if \code{method="mom"},
    and may be omitted for some distributions for which reasonable 
    starting values are computed if \code{method="mle"} (see details).  }
\item{chisqbreaks}{ A numeric vector defining the breaks of the cells used to compute the chi-squared
    statistic. If omitted, these breaks are automatically computed from the data 
    in order to reach roughly the same number of observations per cell, roughly equal to the argument
    \code{meancount}, or sligthly more if there are some ties. }
\item{meancount}{ The mean number of observations per cell expected for the definition of the breaks of the cells 
     used to compute the chi-squared statistic. This argument will not be taken into account if the breaks
     are directly defined in the argument \code{chisqbreaks}. If \code{chisqbreaks} and \code{meancount} are both
     omitted, \code{meancount} is fixed in order to obtain roughly $(4n)^{2/5}$ cells with n the length of the dataset. }
\item{x}{ an object of class 'fitdist'.  }
\item{object}{ an object of class 'fitdist'.  }
\item{breaks}{ If \code{"default"} the histogram is plotted with the function \code{hist} 
    with its default breaks definition. Else \code{breaks} is passed to the function \code{hist}.
    This argument is not taken into account with discrete distributions: \code{"binom"}, 
    \code{"nbinom"}, \code{"geom"}, \code{"hyper"} and \code{"pois"}. }
\item{...}{ further arguments passed to or from other methods }
}

\details{
    When \code{method="mle"},
    maximum likelihood estimations of the distribution parameters are computed using 
    the function \code{\link{mledist}}. 
    Direct optimization of the log-likelihood is performed using \code{\link{optim}}, with its
    default method \code{"Nelder-Mead"} for distributions characterized by more than one parameter and
    the method \code{"BFGS"} for distributions characterized by only one parameter. 
    For the following named distributions, reasonable starting values will 
    be computed if \code{start} is omitted : \code{"norm"}, \code{"lnorm"},
    \code{"exp"} and \code{"pois"}, \code{"cauchy"}, \code{"gamma"}, \code{"logis"},
    \code{"nbinom"} (parametrized by mu and size), \code{"geom"}, \code{"beta"} and \code{"weibull"}. Note that these starting 
    values may not be good enough if the fit is poor. The function is not able to fit a uniform distribution.     
    With the parameter estimates, the function returns the log-likelihood and the standard errors of 
    the estimates calculated from the 
    Hessian at the solution found by \code{optim}.

    
    When \code{method="mom"},
    the estimated values of the distribution parameters are provided only for the following
    distributions : \code{"norm"}, \code{"lnorm"}, \code{"pois"}, \code{"exp"}, \code{"gamma"},
    \code{"nbinom"}, \code{"geom"}, \code{"beta"}, \code{"unif"} and \code{"logis"}.
    For distributions characterized by one parameter (\code{"geom"}, \code{"pois"} and \code{"exp"}), this parameter is simply 
    estimated by matching theoretical and observed means, and for distributions characterized by 
    two parameters, these parameters are estimated by matching theoretical and observed means
    and variances (Vose, 2000). 
     
    Goodness-of-fit statistics are computed. The Chi-squared statistic is computed using cells defined by the argument 
    \code{chisqbreaks} or cells automatically defined from the data in order 
    to reach roughly the same number of observations per cell, roughly equal to the argument
    \code{meancount}, or sligthly more if there are some ties. If \code{chisqbreaks} and \code{meancount} are both
     omitted, \code{meancount} is fixed in order to obtain roughly $(4n)^{2/5}$ cells, with n the length of the dataset (Vose, 2000).
     The Chi-squared statistic is not computed if the program fails
    to define enough cells due to a too small dataset. When the Chi-squared statistic is computed, and if the degree
    of freedom (nb of cells - nb of parameters - 1)  of the corresponding distribution is strictly positive, the p-value
    of the Chi-squared test is returned.
    
    For the distributions assumed continuous (all but \code{"binom"}, 
    \code{"nbinom"}, \code{"geom"}, \code{"hyper"} and \code{"pois"}), Kolmogorov-Smirnov and Anderson-Darling  
    statistics are also computed, as defined by Cullen and Frey (1999). 
    
    An approximate Kolmogorov-Smirnov test is 
    performed by assuming the distribution parameters known. The critical value defined by Stephens (1986) 
    for a completely specified distribution is used to reject or not the 
    distribution at the significance level 0.05. Because of this approximation, the result of the test
    (decision of rejection of the distribution or not) is returned only for datasets with more than 30 observations. 
    Note that this approximate test may be too conservative. 
    
    For datasets with more than 5 observations and for distributions for 
    which the test is described by Stephens (1986) (\code{"norm"}, \code{"lnorm"},
    \code{"exp"}, \code{"cauchy"}, \code{"gamma"}, \code{"logis"} and \code{"weibull"}),
    the Anderson-darling test is performed as described by Stephens (1986). This test takes into account the
    fact that the parameters are not known but estimated from the data. The result is the decision to reject
    or not the distribution at the significance level 0.05.
    
    The plot of an object of class "fitdist" returned by \code{fitdist} uses the function \code{\link{plotdist}}.
}

\value{ 
    \code{fitdist} returns an object of class 'fitdist', a list with 16 components,
    \item{ estimate }{ the parameter estimates }
    \item{ method }{ the character string coding for the fitting method : 
        \code{"mle"} for 'maximum likelihood' and \code{"mom"} for 'matching moments' }
    \item{ sd }{ the estimated standard errors or \code{NULL} if \code{method="mom"} }
    \item{ cor }{ the estimated correlation matrix or \code{NULL} if \code{method="mom"}}
    \item{ loglik }{ the log-likelihood or \code{NULL} if \code{method="mom"}}
    \item{ n }{ the length of the data set }
    \item{ data }{ the dataset }
    \item{ distname }{ the name of the distribution }
    \item{ chisq }{ the Chi-squared statistic or \code{NULL} if not computed }
    \item{ chisqbreaks }{ breaks used to define cells in the Chi-squared statistic }
    \item{ chisqpvalue }{ p-value of the Chi-squared statistic or \code{NULL} if not computed}
    \item{ chisqdf }{ degree of freedom of the Chi-squared distribution or \code{NULL} if not computed }
    \item{ chisqtable }{ a table with observed and theoretical counts used for the Chi-squared calculations }
    \item{ ad }{ the Anderson-Darling statistic or \code{NULL} if not computed }
    \item{ adtest }{ the decision of the Anderson-Darling test or \code{NULL} if not computed }
    \item{ ks }{ the Kolmogorov-Smirnov statistic or \code{NULL} if not computed }
    \item{ kstest }{ the decision of the Kolmogorov-Smirnov test or \code{NULL} if not computed }
}

\seealso{ 
    \code{\link{plotdist}}, \code{\link{optim}}, \code{\link{mledist}}, \code{\link{momdist}} 
    and \code{\link{fitdistcens}}.
}

\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 81-155.

Stephens MA (1986) Tests based on edf statistics. In Goodness-of-fit techniques (D'Agostino RB and
Stephens MA, eds), Marcel dekker, New York, pp. 97-194.

Venables WN and Ripley BD (2002) Modern applied statistics with S. Springer, New York, pp. 435-446.

Vose D (2000) Risk analysis, a quantitative guide. John Wiley & Sons Ltd, Chischester, England, pp. 99-143.
}


\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vet-lyon.fr} 
}

%\note{  }

\examples{
x1<-c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
f1<-fitdist(x1,"norm")
print(f1)
plot(f1)
summary(f1)
f1$chisqtable

f1b<-fitdist(x1,"norm",method="mom",meancount=6)
summary(f1b)
f1b$chisqtable

f1c<-fitdist(x1,"lnorm",method="mom",meancount=6)
summary(f1c)
f1c$chisqtable

dgumbel<-function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel<-function(q,a,b) exp(-exp((a-q)/b))
qgumbel<-function(p,a,b) a-b*log(-log(p))
f1c<-fitdist(x1,"gumbel",start=list(a=10,b=5))
print(f1c)
plot(f1c)

x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
f2<-fitdist(x2,"pois",chisqbreaks=c(0,1))
plot(f2)
summary(f2)
f2$chisqtable

xw<-rweibull(n=100,shape=2,scale=1)
fa<-fitdist(xw,"weibull")
summary(fa)
fa$chisqtable
fb<-fitdist(xw,"gamma")
summary(fb)
fc<-fitdist(xw,"exp")
summary(fc)


}

\keyword{ distribution }% at least one, from doc/KEYWORDS
