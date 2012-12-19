\name{fitdist}
\alias{fitdist}
\alias{plot.fitdist}
\alias{print.fitdist}
\alias{summary.fitdist}
\alias{fitdistrplus}
\title{ Fit of univariate distributions to non-censored data}

\description{
 Fit of univariate distributions to non-censored data by maximum likelihood, 
 quantile matching or moment matching. Generic methods are \code{print}, \code{plot}, 
 \code{summary} and \code{quantile}.
}

\usage{
fitdist(data, distr, method = c("mle", "mme", "qme", "mge"), 
    start=NULL, fix.arg=NULL, \dots)
    
\method{print}{fitdist}(x, \dots)

\method{plot}{fitdist}(x, breaks="default", \dots)

\method{summary}{fitdist}(object, \dots)

}

\arguments{
\item{data}{A numeric vector.}
\item{distr}{A character string \code{"name"} naming a distribution for which the corresponding
    density function \code{dname}, the corresponding distribution function \code{pname} and the 
    corresponding quantile function \code{qname} must be defined, or directly the density function.}
\item{method}{A character string coding for the fitting method: 
    \code{"mle"} for 'maximum likelihood estimation', \code{"mme"} for 'moment matching estimation', 
    \code{"qme"} for 'quantile matching estimation' and \code{"mge"} for 'maximum goodness-of-fit estimation'.}
\item{start}{A named list giving the initial values of parameters of the named distribution. 
    This argument may be omitted for some distributions for which reasonable 
    starting values are computed (see details), 
    and will not be taken into account if a closed formula is used to estimate parameters.  }
\item{fix.arg}{An optional named list giving the values of parameters of the named distribution
    that must kept fixed rather than estimated. 
    The use of this argument is not possible if \code{method="mme"} and a closed formula is used. }
\item{x}{An object of class 'fitdist'.  }
\item{object}{An object of class 'fitdist'.  }
\item{breaks}{If \code{"default"} the histogram is plotted with the function \code{hist} 
    with its default breaks definition. Else \code{breaks} is passed to the function \code{hist}.
    This argument is not taken into account with discrete distributions: \code{"binom"}, 
    \code{"nbinom"}, \code{"geom"}, \code{"hyper"} and \code{"pois"}. }
\item{\dots}{Further arguments to be passed to generic functions, or to one of the functions
    \code{"mledist"},  \code{"mmedist"}, \code{"qmedist"} or \code{"mgedist"}
    depending of the chosen method
    (see the help pages of these functions for details).}
}

\details{
    Fittting methods:
\describe{
    \item{When \code{method="mle"}}{
    maximum likelihood estimations of the distribution parameters are computed using 
    the function \code{\link{mledist}}. 
    }
    \item{When \code{method="mme"}}{
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
    }
    \item{When \code{method = "qme"}}{
    the function carries out the quantile matching numerically, by minimization of the
    sum of squared differences between observed and theoretical quantiles.
    The use of this method requires an additional argument \code{probs}, 
    defined as the numeric vector of the probabilities 
    for which the quantile matching is done, of length equal to the number of parameters to estimate
    (see \code{\link{qmedist}} for details).
    }
    \item{When \code{method = "mge"}}{
    the distribution parameters are estimated by maximization  
    of goodness-of-fit (or minimization of a goodness-of-fit distance). The use of this method 
    requires an additional argument 
    \code{gof} coding for the goodness-of-fit distance chosen. 
    One may use the classical Cramer-von Mises distance (\code{"CvM"}), the classical
    Kolmogorov-Smirnov distance (\code{"KS"}), the classical Anderson-Darling distance (\code{"AD"})
    which gives more weight to the tails of the distribution, 
    or one of the variants of this last distance proposed by Luceno (2006)
    (see \code{\link{mgedist}} for more details). This method is not suitable for discrete distributions.
    }
}
        
    By default direct optimization of the log-likelihood (or other criteria depending
    of the chosen method) is performed using \code{\link{optim}}, 
    with the "Nelder-Mead" method for distributions characterized by more than one parameter
    and the "BFGS" method for distributions characterized by only one parameter. 
    The method used in \code{\link{optim}} may be chosen or another optimization method
    may be chosen using \dots argument (see \code{\link{mledist}} for details).
    For the following named distributions, reasonable starting values will 
    be computed if \code{start} is omitted : \code{"norm"}, \code{"lnorm"}, 
    \code{"exp"} and \code{"pois"}, \code{"cauchy"}, \code{"gamma"}, \code{"logis"}, 
    \code{"nbinom"} (parametrized by mu and size), \code{"geom"}, \code{"beta"}, 
    \code{"unif"}, and \code{"weibull"}. 
    Note that these starting 
    values may not be good enough if the fit is poor.      
    With the parameter estimates, the function returns the log-likelihood whatever the estimation method
    and for maximum likelihood estimation the standard errors of 
    the estimates calculated from the 
    Hessian at the solution found by \code{optim} or by the user-supplied function passed to mledist.

    NB: if your data values are particularly small or large, a scaling may be needed 
    before the optimization process. See example (14).
    
}

\value{ 
    \code{fitdist} returns an object of class 'fitdist', a list with the following components:
    \item{estimate }{ the parameter estimates }
    \item{method }{ the character string coding for the fitting method : 
        \code{"mle"} for 'maximum likelihood estimation', \code{"mme"} for 'matching moment estimation',
         \code{"qme"} for 'matching quantile estimation'
         and \code{"mge"} for 'maximum goodness-of-fit estimation' }
    \item{sd}{ the estimated standard errors or \code{NULL} if not available }
    \item{cor}{ the estimated correlation matrix or \code{NULL} if not available}
    \item{loglik}{ the log-likelihood}
    \item{aic}{ the Akaike information criterion}
    \item{bic}{ the the so-called BIC or SBC (Schwarz Bayesian criterion)}
    \item{n}{ the length of the data set }
    \item{data}{ the dataset }
    \item{distname}{ the name of the distribution }
    \item{fix.arg}{ the named list giving the values of parameters of the named distribution
        that must kept fixed rather than estimated by maximum likelihood or \code{NULL} 
        if there are no such parameters. }
    \item{dots}{ the list of  further arguments passed in \dots to be used in \code{bootdist} 
        in iterative calls to \code{mledist}, \code{mmedist}, \code{qmedist}, \code{mgedist} or 
        \code{NULL} if no such arguments.}
        

    Generic functions:
\describe{       
    \item{\code{print}}{
    The print of a 'fitdist' object shows few traces about the fitting method and the fitted distribution.
    }
    \item{\code{summary}}{
    The summary provides the parameter estimates of the fitted distribution, the log-likelihood, 
    AIC and BIC statistics and when the maximum likelihood is used, the standard errors of the 
    parameter estimates and the correlation matrix between parameter estimates.
    }
    \item{\code{plot}}{
    The plot of an object of class "fitdist" returned by \code{fitdist} uses the function 
    \code{\link{plotdist}}. An object of class "fitdist" or a list of objects of class
    "fitdist" corresponding to various fits using the same data set may also be plotted
    using a cdf plot (function \code{\link{cdfcomp}}),
    a density plot(function \code{\link{denscomp}}),
    a density Q-Q plot (function \code{\link{qqcomp}}),
    or a P-P plot (function \code{\link{ppcomp}}).
    }
    
}

}

\seealso{ 
    See \code{\link{mledist}}, \code{\link{mmedist}}, \code{\link{qmedist}}, 
    \code{\link{mgedist}} for details on parameter estimation.
    See \code{\link{gofstat}} for goodness-of-fit statistics.
    See \code{\link{plotdist}}, \code{\link{graphcomp}} for graphs.
    See \code{\link{bootdist}} for bootstrap procedures 
    and \code{\link{fitdistcens}} for censored-data fitting methods.
    See \code{\link{optim}} for base R optimization procedures.
    See \code{\link{quantile.fitdist}} for another generic function to calculate 
    quantiles from the fitted distribution.
}

\references{ 
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure assessment. Plenum Press, USA, pp. 81-155.

Venables WN and Ripley BD (2002) Modern applied statistics with S. Springer, New York, pp. 435-446.

Vose D (2000) Risk analysis, a quantitative guide. John Wiley & Sons Ltd, Chischester, England, pp. 99-143.
}


\author{ 
Marie-Laure Delignette-Muller \email{marielaure.delignettemuller@vetagro-sup.fr} 
and Christophe Dutang.
}

\examples{

# (1) fit of a gamma distribution by maximum likelihood estimation
#

data(groundbeef)
serving <- groundbeef$serving
fitg <- fitdist(serving, "gamma")
summary(fitg)
plot(fitg)
cdfcomp(fitg, addlegend=FALSE)
denscomp(fitg, addlegend=FALSE)
ppcomp(fitg, addlegend=FALSE)
qqcomp(fitg, addlegend=FALSE)


# (2) use the moment matching estimation (using a closed formula)
#

fitgmme <- fitdist(serving, "gamma", method="mme")
summary(fitgmme)

# (3) Comparison of various fits 
#

fitW <- fitdist(serving, "weibull")
fitg <- fitdist(serving, "gamma")
fitln <- fitdist(serving, "lnorm")
summary(fitW)
summary(fitg)
summary(fitln)
cdfcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))
denscomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))
qqcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))
ppcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))

# (4) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view 
# dedicated to probability distributions
#

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

fitgumbel <- fitdist(serving, "gumbel", start=list(a=10, b=10))
summary(fitgumbel)
plot(fitgumbel)

# (5) fit discrete distributions (Poisson and negative binomial)
#

data(toxocara)
number <- toxocara$number
fitp <- fitdist(number,"pois")
summary(fitp)
plot(fitp)
fitnb <- fitdist(number,"nbinom")
summary(fitnb)
plot(fitnb)

# (6) how to change the optimisation method?
#

data(groundbeef)
serving <- groundbeef$serving
fitdist(serving, "gamma", optim.method="Nelder-Mead")
fitdist(serving, "gamma", optim.method="BFGS") 
fitdist(serving, "gamma", optim.method="SANN")

# (7) custom optimization function
#

#create the sample
set.seed(1234)
mysample <- rexp(100, 5)
mystart <- list(rate=8)

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


# (8) custom optimization function - another example with the genetic algorithm
#
\dontrun{
    #set a sample
    fit1 <- fitdist(serving, "gamma")
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
    fit2 <- fitdist(serving, "gamma", custom.optim=mygenoud, nvars=2,    
        Domains=cbind(c(0, 0), c(10, 10)), boundary.enforcement=1, 
        print.level=1, hessian=TRUE)

    summary(fit2)
}

# (9) estimation of the standard deviation of a gamma distribution 
# by maximum likelihood with the shape fixed at 4 using the argument fix.arg
#

data(groundbeef)
serving <- groundbeef$serving
f1c  <- fitdist(serving,"gamma",start=list(rate=0.1),fix.arg=list(shape=4))
summary(f1c)
plot(f1c)

# (10) fit of a Weibull distribution to serving size data 
# by maximum likelihood estimation
# or by quantile matching estimation (in this example 
# matching first and third quartiles)
#

data(groundbeef)
serving <- groundbeef$serving
fWmle <- fitdist(serving, "weibull")
summary(fWmle)
plot(fWmle)
gofstat(fWmle)

fWqme <- fitdist(serving, "weibull", method="qme", probs=c(0.25, 0.75))
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
    fP <- fitdist(x4, "pareto", method="mme", order=c(1, 2), memp="memp", 
    start=c(shape=10, scale=10), lower=1, upper=Inf)
    summary(fP)
    plot(fP)

}

# (12) Fit of a Weibull distribution to serving size data by maximum 
# goodness-of-fit estimation using all the distances available
# 

data(groundbeef)
serving <- groundbeef$serving
(f1 <- fitdist(serving, "weibull", method="mge", gof="CvM"))
(f2 <- fitdist(serving, "weibull", method="mge", gof="KS"))
(f3 <- fitdist(serving, "weibull", method="mge", gof="AD"))
(f4 <- fitdist(serving, "weibull", method="mge", gof="ADR"))
(f5 <- fitdist(serving, "weibull", method="mge", gof="ADL"))
(f6 <- fitdist(serving, "weibull", method="mge", gof="AD2R"))
(f7 <- fitdist(serving, "weibull", method="mge", gof="AD2L"))
(f8 <- fitdist(serving, "weibull", method="mge", gof="AD2"))
cdfcomp(list(f1, f2, f3, f4, f5, f6, f7, f8))
cdfcomp(list(f1, f2, f3, f4, f5, f6, f7, f8), xlogscale=TRUE, xlim=c(8, 250), verticals=TRUE)
denscomp(list(f1, f2, f3, f4, f5, f6, f7, f8))

# (13) Fit of a uniform distribution using maximum likelihood 
# (a closed formula is used in this special case where the loglikelihood is not defined),
# or maximum goodness-of-fit with Cramer-von Mises or Kolmogorov-Smirnov distance
# 

set.seed(1234)
u <- runif(50, min=5, max=10)

fumle <- fitdist(u, "unif", method="mle")
summary(fumle)
plot(fumle)
gofstat(fumle)

fuCvM <- fitdist(u, "unif", method="mge", gof="CvM")
summary(fuCvM)
plot(fuCvM)
gofstat(fuCvM)

fuKS <- fitdist(u, "unif", method="mge", gof="KS")
summary(fuKS)
plot(fuKS)
gofstat(fuKS)

# (14) scaling problem
# the simulated dataset (below) has particularly small values, hence without scaling (10^0),
# the optimization raises an error. The for loop shows how scaling by 10^i
# for i=1,...,6 makes the fitting procedure work correctly.

set.seed(1234)
x2 <- rnorm(100, 1e-4, 2e-4)

for(i in 0:6)
        cat(i, try(fitdist(x2*10^i, "cauchy", method="mle")$estimate, silent=TRUE), "\n")

# (15) Fit of a lognormal distribution on acute toxicity values of endosulfan for
# nonarthropod invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5 percent quantile value of 
# the fitted distribution, what is called the 5 percent hazardous concentration (HC5)
# in ecotoxicology, with its two-sided 95 percent confidence interval calculated by 
# parametric bootstrap or non parametric bootstrap, and
# its one-sided (greater) 95 percent confidence interval calculated by 
# parametric bootstrap
#
data(endosulfan)
ATV <- subset(endosulfan, group == "NonArthroInvert")$ATV
log10ATV <- log10(subset(endosulfan, group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV, "norm")

quantile(fln, probs = 0.05)
quantile(fln, probs = c(0.05, 0.1, 0.2))


# (16) Fit of a triangular distribution using Cramer-von Mises or
# Kolmogorov-Smirnov distance
# 

\dontrun{
set.seed(1234)
require(mc2d)
t <- rtriang(100, min=5, mode=6, max=10)
fCvM <- fitdist(t, "triang", method="mge", start = list(min=4, mode=6,max=9), gof="CvM")
fKS <- fitdist(t, "triang", method="mge", start = list(min=4, mode=6,max=9), gof="KS")
cdfcomp(list(fCvM,fKS))
}

}
\keyword{ distribution }
