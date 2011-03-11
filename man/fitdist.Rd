\name{fitdist}
\alias{fitdist}
\alias{plot.fitdist}
\alias{print.fitdist}
\alias{summary.fitdist}
\alias{fitdistrplus}
\title{ Fit of univariate distributions to non-censored data}

\description{
 Fit of univariate distributions to non-censored data by maximum likelihood,
 quantile matching or moment matching.  
}

\usage{
fitdist(data, distr, method=c("mle", "mme", "qme", "mge"), 
    start=NULL, fix.arg=NULL,  ...) 
\method{print}{fitdist}(x,...)
\method{plot}{fitdist}(x,breaks="default",...)
\method{summary}{fitdist}(object,...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
\item{data}{ A numeric vector.  }
\item{distr}{ A character string \code{"name"} naming a distribution for which the corresponding
    density function \code{dname}, the corresponding distribution function \code{pname} and the 
    corresponding quantile function \code{qname} must be defined, or directly the density function.}
\item{method}{ A character string coding for the fitting method: 
    \code{"mle"} for 'maximum likelihood estimation', \code{"mme"} for 'moment matching estimation',
    \code{"qme"} for 'quantile matching estimation' and \code{"mge"} for 'maximum goodness-of-fit estimation'.}
\item{start}{ An named list giving the initial values of parameters of the named distribution. 
    This argument may be omitted for some distributions for which reasonable 
    starting values are computed (see details), 
    and will not be taken into account if a closed formula is used to estimate parameters.  }
\item{fix.arg}{ An optional named list giving the values of parameters of the named distribution
    that must kept fixed rather than estimated. 
    The use of this argument is not possible if \code{method="mme"} and a closed formula is used. }
\item{x}{ an object of class 'fitdist'.  }
\item{object}{ an object of class 'fitdist'.  }
\item{breaks}{ If \code{"default"} the histogram is plotted with the function \code{hist} 
    with its default breaks definition. Else \code{breaks} is passed to the function \code{hist}.
    This argument is not taken into account with discrete distributions: \code{"binom"}, 
    \code{"nbinom"}, \code{"geom"}, \code{"hyper"} and \code{"pois"}. }
\item{\dots}{ further arguments to be passed to generic functions, or to one of the functions
    \code{"mledist"},  \code{"mmedist"}, \code{"qmedist"} or \code{"mgedist"}
    depending of the chosen method
    (see the help pages of these functions for details).}
}

\details{
    When \code{method="mle"},
    maximum likelihood estimations of the distribution parameters are computed using 
    the function \code{\link{mledist}}. 
    
    When \code{method="mme"},
    the estimated values of the distribution parameters are computed by a closed 
    formula for the following distributions : \code{"norm"}, \code{"lnorm"}, 
    \code{"pois"}, \code{"exp"}, \code{"gamma"},
    \code{"nbinom"}, \code{"geom"}, \code{"beta"}, \code{"unif"} and \code{"logis"}.
    For distributions characterized by one parameter (\code{"geom"}, \code{"pois"} and \code{"exp"}), 
    this parameter is simply estimated by matching theoretical and 
    observed means, and for distributions characterized by 
    two parameters, these parameters are estimated by matching theoretical and observed means
    and variances (Vose, 2000). 
    For other distributions, the theoretical and the empirical moments are matched numerically,
    by minimization of the
    sum of squared differences between observed and theoretical moments. In this last case,
    further arguments are needed in the call to \code{fitdist}: \code{order} and \code{memp}
    (see \code{\link{mmedist}} for details).
    
    When \code{method = "qme"},
    the function carries out the quantile matching numerically, by minimization of the
    sum of squared differences between observed and theoretical quantiles.
    The use of this method requires an additional argument \code{probs},
    defined as the numeric vector of the probabilities 
    for which the quantile matching is done, of length equal to the number of parameters to estimate
    (see \code{\link{qmedist}} for details).
    
    When \code{method = "mge"},
    the distribution parameters are estimated by maximization of 
    the goodness-of-fit distance. The use of this method requires an additional argument 
    \code{gof} coding for the goodness-of-fit distance chosen. 
    One may use the classical Cramer-von Mises distance (\code{"CvM"}), the classical
    Kolmogorov-Smirnov distance (\code{"KS"}), the classical Anderson-Darling distance (\code{"AD"})
    which gives more weight to the tails of the distribution,
    or one of the variants of this last distance proposed by Luceno (2006)
    (see \code{\link{mgedist}} for more details). This method is not suitable for discrete distributions.

    By default direct minimization of the log-likelihood (or other criteria depending
    of the chosen method) is performed using \code{\link{optim}},
    with the "Nelder-Mead" method for distributions characterized by more than one parameter
    and the "BFGS" method for distributions characterized by only one parameter. 
    The method used in \code{\link{optim}} may be chosen or another optimization method
    may be chosen using ... argument (see \code{\link{mledist}} for details).
    For the following named distributions, reasonable starting values will 
    be computed if \code{start} is omitted : \code{"norm"}, \code{"lnorm"},
    \code{"exp"} and \code{"pois"}, \code{"cauchy"}, \code{"gamma"}, \code{"logis"},
    \code{"nbinom"} (parametrized by mu and size), \code{"geom"}, \code{"beta"} and \code{"weibull"}. 
    Note that these starting 
    values may not be good enough if the fit is poor. The function is not able to fit a uniform distribution.     
    With the parameter estimates, the function returns the log-likelihood whatever the estimation method
    and for maximum likelihood estimation the standard errors of 
    the estimates calculated from the 
    Hessian at the solution found by \code{optim} or by the user-supplied function passed to mledist.

         
    The plot of an object of class "fitdist" returned by \code{fitdist} uses the function 
    \code{\link{plotdist}}.
}

\value{ 
    \code{fitdist} returns an object of class 'fitdist', a list with following components,
    \item{ estimate }{ the parameter estimates }
    \item{ method }{ the character string coding for the fitting method : 
        \code{"mle"} for 'maximum likelihood estimation', \code{"mme"} for 'matching moment estimation'
        and \code{"qme"} for 'matching quantile estimation' }
    \item{ sd }{ the estimated standard errors or \code{NULL} if not available }
    \item{ cor }{ the estimated correlation matrix or \code{NULL} if not available}
    \item{ loglik }{ the log-likelihood}
     \item{ aic }{ the Akaike information criterion}
    \item{ bic }{ the the so-called BIC or SBC (Schwarz Bayesian criterion)}
   \item{ n }{ the length of the data set }
    \item{ data }{ the dataset }
    \item{ distname }{ the name of the distribution }
    \item{ fix.arg }{ the named list giving the values of parameters of the named distribution
    that must kept fixed rather than estimated by maximum likelihood or NULL if there are no such parameters. }
    \item{ dots }{ the list of  further arguments passed in \dots to be used in \code{bootdist} 
    in iterative calls to \code{mledist},
    \code{mmedist}, \code{qmedist}, \code{mgedist} or NULL if no such arguments}
}

\seealso{ 
    \code{\link{plotdist}}, \code{\link{optim}}, \code{\link{mledist}}, \code{\link{mmedist}}, \code{\link{qmedist}},
    \code{\link{mgedist}}, \code{\link{gofstat}}
    and \code{\link{fitdistcens}}.
}

\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 81-155.

Venables WN and Ripley BD (2002) Modern applied statistics with S. Springer, New York, pp. 435-446.

Vose D (2000) Risk analysis, a quantitative guide. John Wiley & Sons Ltd, Chischester, England, pp. 99-143.
}


\author{ 
Marie-Laure Delignette-Muller \email{ml.delignette@vetagro-sup.fr} and
Christophe Dutang
}

%\note{  }

\examples{

# (1) basic fit of a normal distribution with maximum likelihood estimation
#

x1 <- c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,
13.2,8.4,6.3,8.9,5.2,10.9,14.4)
f1 <- fitdist(x1,"norm")
print(f1)
plot(f1)
summary(f1)
gofstat(f1)

# (2) use the moment matching estimation (using a closed formula)
#

f1b <- fitdist(x1,"norm",method="mme")
summary(f1b)

# (3) moment matching estimation (using a closed formula) 
# for log normal distribution
#

f1c <- fitdist(x1,"lnorm",method="mme")
summary(f1c)

# (4) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view 
# dedicated to probability distributions
#

dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q,a,b) exp(-exp((a-q)/b))
qgumbel <- function(p,a,b) a-b*log(-log(p))

f1c <- fitdist(x1,"gumbel",start=list(a=10,b=5))
print(f1c)
plot(f1c)

# (5) fit a discrete distribution (Poisson)
#

x2<-c(rep(4,1),rep(2,3),rep(1,7),rep(0,12))
f2<-fitdist(x2,"pois")
plot(f2)
summary(f2)
gofstat(f2)

# (6) how to change the optimisation method?
#

fitdist(x1,"gamma",optim.method="Nelder-Mead")
fitdist(x1,"gamma",optim.method="BFGS") 
fitdist(x1,"gamma",optim.method="L-BFGS-B",lower=c(0,0))
fitdist(x1,"gamma",optim.method="SANN")

# (7) custom optimization function
#

#create the sample
mysample <- rexp(100, 5)
mystart <- 8

res1 <- fitdist(mysample, dexp, start= mystart, optim.method="Nelder-Mead")

#show the result
summary(res1)

#the warning tell us to use optimise, because the Nelder-Mead is not adequate.

#to meet the standard 'fn' argument and specific name arguments, we wrap optimize,
myoptimize <- function(fn, par, ...) 
{
    res <- optimize(f=fn, ..., maximum=FALSE)  
    #assume the optimization function minimize
    
    standardres <- c(res, convergence=0, value=res$objective, 
        par=res$minimum, hessian=NA)
    
    return(standardres)
}

#call fitdist with a 'custom' optimization function
res2 <- fitdist(mysample, dexp, start=mystart, custom.optim=myoptimize, 
    interval=c(0, 100))

#show the result
summary(res2)


# (8) custom optimisation function - another example with the genetic algorithm
#
\dontrun{
    #set a sample
    x1 <- c(6.4, 13.3, 4.1, 1.3, 14.1, 10.6, 9.9, 9.6, 15.3, 22.1,
         13.4, 13.2, 8.4, 6.3, 8.9, 5.2, 10.9, 14.4) 
    fit1 <- fitdist(x1, "gamma")
    summary(fit1)

    #wrap genoud function rgenoud package
    mygenoud <- function(fn, par, ...) 
    {
        require(rgenoud)
        res <- genoud(fn, starting.values=par, ...)        
        standardres <- c(res, convergence=0)
            
        return(standardres)
    }

    #call fitdist with a 'custom' optimization function
    fit2 <- fitdist(x1, "gamma", custom.optim=mygenoud, nvars=2,    
        Domains=cbind(c(0,0), c(10, 10)), boundary.enforcement=1, 
        print.level=1, hessian=TRUE)

    summary(fit2)
}

# (9) estimation of the standard deviation of a normal distribution 
# by maximum likelihood with the mean fixed at 10 using the argument fix.arg
#
fitdist(x1,"norm",start=list(sd=5),fix.arg=list(mean=10))

# (10) fit of a Weibull distribution to serving size data by maximum likelihood estimation
#  or by quantile matching estimation (in this example matching first and third quartiles)
#
data(groundbeef)
serving <- groundbeef$serving

fWmle <- fitdist(serving,"weibull")
summary(fWmle)
plot(fWmle)
gofstat(fWmle)

fWqme <- fitdist(serving,"weibull",method="qme",probs=c(0.25,0.75))
summary(fWqme)
plot(fWqme)
gofstat(fWqme)


# (11) Fit of a Pareto distribution by numerical moment matching estimation
#
\dontrun{
    require(actuar)
    #simulate a sample
    x4 <- rpareto(1000, 6, 2)

    #empirical raw moment
    memp <- function(x, order)
        ifelse(order == 1, mean(x), sum(x^order)/length(x))


    #fit
    fP <- fitdist(x4, "pareto", method="mme",order=c(1, 2), memp="memp", 
    start=c(10, 10), lower=1, upper=Inf)
    summary(fP)

}

# (12) Fit of a Weibull distribution to serving size data by maximum 
# goodness-of-fit estimation using all the distances available
# 

data(groundbeef)
serving <- groundbeef$serving
fitdist(serving,"weibull",method="mge",gof="CvM")
fitdist(serving,"weibull",method="mge",gof="KS")
fitdist(serving,"weibull",method="mge",gof="AD")
fitdist(serving,"weibull",method="mge",gof="ADR")
fitdist(serving,"weibull",method="mge",gof="ADL")
fitdist(serving,"weibull",method="mge",gof="AD2R")
fitdist(serving,"weibull",method="mge",gof="AD2L")
fitdist(serving,"weibull",method="mge",gof="AD2")

# (13) Fit of a uniform distribution using Cramer-von Mises or
# Kolmogorov-Smirnov distance
# 

u <- runif(50,min=5,max=10)

fuCvM <- fitdist(u,"unif",method="mge",gof="CvM")
summary(fuCvM)
plot(fuCvM)
gofstat(fuCvM)

fuKS <- fitdist(u,"unif",method="mge",gof="KS")
summary(fuKS)
plot(fuKS)
gofstat(fuKS)


}
\keyword{ distribution }
