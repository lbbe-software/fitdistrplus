\name{fitdist}
\alias{fitdist}
\alias{plot.fitdist}
\alias{print.fitdist}
\alias{summary.fitdist}
\alias{logLik.fitdist}
\alias{vcov.fitdist}
\alias{coef.fitdist}
\title{ Fit of univariate distributions to non-censored data}

\description{
 Fit of univariate distributions to non-censored data by maximum likelihood (mle), 
 moment matching (mme), quantile matching (qme) or 
 maximizing goodness-of-fit estimation (mge).
 The latter is also known as minimizing distance estimation.
 Generic methods are \code{print}, \code{plot}, 
\code{summary}, \code{quantile}, \code{logLik}, \code{vcov} and \code{coef}.
}

\usage{
fitdist(data, distr, method = c("mle", "mme", "qme", "mge"), 
    start=NULL, fix.arg=NULL, discrete, keepdata = TRUE, keepdata.nb=100, \dots)
    
\method{print}{fitdist}(x, \dots)

\method{plot}{fitdist}(x, breaks="default", \dots)

\method{summary}{fitdist}(object, \dots)

\method{logLik}{fitdist}(object, \dots)

\method{vcov}{fitdist}(object, \dots)

\method{coef}{fitdist}(object, \dots)


}

\arguments{
\item{data}{A numeric vector.}
\item{distr}{A character string \code{"name"} naming a distribution for which the corresponding
    density function \code{dname}, the corresponding distribution function \code{pname} and the 
    corresponding quantile function \code{qname} must be defined, or directly the density function.}
\item{method}{A character string coding for the fitting method: 
    \code{"mle"} for 'maximum likelihood estimation', \code{"mme"} for 'moment matching estimation', 
    \code{"qme"} for 'quantile matching estimation' and \code{"mge"} for 'maximum goodness-of-fit estimation'.}
\item{start}{A named list giving the initial values of parameters of the named distribution
    or a function of data computing initial values and returning a named list.
    This argument may be omitted (default) for some distributions for which reasonable 
    starting values are computed (see the 'details' section of  \code{\link{mledist}}). 
    It may not be into account for closed-form formulas.}
\item{fix.arg}{An optional named list giving the values of fixed parameters of the named distribution
    or a function of data computing (fixed) parameter values and returning a named list. 
    Parameters with fixed value are thus NOT estimated by this maximum likelihood procedure.
    The use of this argument is not possible if \code{method="mme"} and a closed-form formula is used.}
\item{x}{An object of class \code{"fitdist"}.  }
\item{object}{An object of class \code{"fitdist"}.  }
\item{breaks}{If \code{"default"} the histogram is plotted with the function \code{hist} 
    with its default breaks definition. Else \code{breaks} is passed to the function \code{hist}.
    This argument is not taken into account with discrete distributions: \code{"binom"}, 
    \code{"nbinom"}, \code{"geom"}, \code{"hyper"} and \code{"pois"}. }
    
\item{keepdata}{a logical. If \code{TRUE}, dataset is returned, 
    otherwise only a sample subset is returned.}
\item{keepdata.nb}{When \code{keepdata=FALSE}, the length (>1) of the subset returned.}
\item{discrete}{ If TRUE, the distribution is considered as discrete.
    If \code{discrete} is missing, 
	  \code{discrete} is automaticaly set to \code{TRUE} when \code{distr} belongs to 
	  \code{"binom"}, \code{"nbinom"}, \code{"geom"},
    \code{"hyper"} or \code{"pois"} and to \code{FALSE} in the other cases. It is thus recommended
    to enter this argument when using another discrete distribution. This argument will not directly affect
    the results of the fit but will be passed to functions \code{\link{gofstat}}, \code{\link{plotdist}} 
    and \code{\link{cdfcomp}}.  }
\item{\dots}{Further arguments to be passed to generic functions, or to one of the functions
    \code{"mledist"},  \code{"mmedist"}, \code{"qmedist"} or \code{"mgedist"}
    depending of the chosen method. See \code{\link{mledist}}, \code{\link{mmedist}}, \code{\link{qmedist}}, 
    \code{\link{mgedist}} for details on parameter estimation.}
}

\details{

    It is assumed that the \code{distr} argument specifies the distribution by the
    probability density function, the cumulative distribution function and
    the quantile function (d, p, q).

    The four possible fitting methods are described below:
\describe{
    \item{When \code{method="mle"}}{
    Maximum likelihood estimation consists in maximizing the log-likelihood.
    A numerical optimization is carried out in \code{\link{mledist}} via \code{optim}
    to find the best values (see \code{\link{mledist}} for details).
    }
    \item{When \code{method="mme"}}{
    Moment matching estimation consists in equalizing theoretical and empirical moments.
    Estimated values of the distribution parameters are computed by a closed-form 
    formula for the following distributions : \code{"norm"}, \code{"lnorm"}, 
    \code{"pois"}, \code{"exp"}, \code{"gamma"}, 
    \code{"nbinom"}, \code{"geom"}, \code{"beta"}, \code{"unif"} and \code{"logis"}.
    Otherwise the theoretical and the empirical moments are matched numerically, 
    by minimization of the sum of squared differences between observed and theoretical moments. 
    In this last case, further arguments are needed in the call to \code{fitdist}: \code{order} and \code{memp}
    (see \code{\link{mmedist}} for details).
    }
    \item{When \code{method = "qme"}}{
    Quantile matching estimation consists in equalizing theoretical and empirical quantile.
    A numerical optimization is carried out in \code{\link{qmedist}} via \code{optim}
    to minimize of the sum of squared differences between observed and theoretical quantiles.
    The use of this method requires an additional argument \code{probs}, 
    defined as the numeric vector of the probabilities for which the quantile(s) is(are) to be matched
    (see \code{\link{qmedist}} for details).
    }
    \item{When \code{method = "mge"}}{
    Maximum goodness-of-fit estimation consists in maximizing a goodness-of-fit statistics.
    A numerical optimization is carried out in \code{\link{mgedist}} via \code{optim}
    to minimize the goodness-of-fit distance. The use of this method 
    requires an additional argument \code{gof} coding for the goodness-of-fit distance chosen. 
    One can use the classical Cramer-von Mises distance (\code{"CvM"}), the classical
    Kolmogorov-Smirnov distance (\code{"KS"}), the classical Anderson-Darling distance (\code{"AD"})
    which gives more weight to the tails of the distribution, 
    or one of the variants of this last distance proposed by Luceno (2006)
    (see \code{\link{mgedist}} for more details). This method is not suitable for discrete distributions.
    }
}
        
    By default, direct optimization of the log-likelihood (or other criteria depending
    of the chosen method) is performed using \code{\link{optim}}, 
    with the "Nelder-Mead" method for distributions characterized by more than one parameter
    and the "BFGS" method for distributions characterized by only one parameter. 
    The optimization algorithm used in \code{\link{optim}} can be chosen or another optimization function
    can be specified using \dots argument (see \code{\link{mledist}} for details).
    \code{start} may be omitted (i.e. \code{NULL}) for some classic distributions
    (see the 'details' section of  \code{\link{mledist}}).
    Note that when errors are raised by \code{optim}, it's a good idea to start by adding traces during
    the optimization process by adding \code{control=list(trace=1, REPORT=1)} in \dots argument.
    
    Once the parameter(s) is(are) estimated, \code{fitdist} computes the log-likelihood 
    for every estimation method and for maximum likelihood estimation the standard errors of 
    the estimates calculated from the Hessian at the solution found by \code{optim} 
    or by the user-supplied function passed to mledist.

    By default (\code{keepdata = TRUE}), the object returned by \code{fitdist} contains 
    the data vector given in input.
    When dealing with large datasets, we can remove the original dataset from the output by
    setting \code{keepdata = FALSE}. In such a case, only \code{keepdata.nb} points (at most) 
    are kept by random subsampling \code{keepdata.nb}-2 points from the dataset and 
    adding the minimum and the maximum. If combined with \code{\link{bootdist}}, and use with
    non-parametric bootstrap be aware that bootstrap is performed on the subset
    randomly selected in \code{fitdist}. Currently, the graphical comparisons of multiple fits 
    is not available in this framework.
    
    Weighted version of the estimation process is available for \code{method = "mle", "mme", "qme"}
    by using \code{weights=\dots}. See the corresponding man page for details.
    Weighted maximum GOF estimation (when \code{method = "mge"}) is not allowed. 
    It is not yet possible to take into account weighths in functions \code{plotdist}, 
    \code{plot.fitdist}, \code{cdfcomp}, \code{denscomp}, \code{ppcomp}, 
    \code{qqcomp}, \code{gofstat} and \code{descdist} (developments planned in the future).

    NB: if your data values are particularly small or large, a scaling may be needed 
    before the optimization process. See example (14) in this man page and 
    examples (14,15) in the test file of the package.
    Please also take a look at the \code{Rmpfr} package available on CRAN for numerical
    accuracy issues.
    
    
}

\value{ 
    \code{fitdist} returns an object of class \code{"fitdist"}, a list with the following components:
    \item{estimate}{ the parameter estimates.}
    \item{method}{ the character string coding for the fitting method : 
        \code{"mle"} for 'maximum likelihood estimation', \code{"mme"} for 'matching moment estimation',
         \code{"qme"} for 'matching quantile estimation'
    and \code{"mge"} for 'maximum goodness-of-fit estimation'.}
    \item{sd}{ the estimated standard errors, \code{NA} if numerically not computable
        or \code{NULL} if not available.}
    \item{cor}{ the estimated correlation matrix, \code{NA} if numerically not computable
        or \code{NULL} if not available.}
    \item{vcov}{ the estimated variance-covariance matrix, \code{NULL} if not available.}
    \item{loglik}{ the log-likelihood.}
    \item{aic}{ the Akaike information criterion.}
    \item{bic}{ the the so-called BIC or SBC (Schwarz Bayesian criterion).}
    \item{n}{ the length of the data set.}
    \item{data}{ the data set.}
    \item{distname}{ the name of the distribution.}
    \item{fix.arg}{ the named list giving the values of parameters of the named distribution
        that must be kept fixed rather than estimated by maximum likelihood or \code{NULL} 
        if there are no such parameters. }
    \item{fix.arg.fun}{the function used to set the value of \code{fix.arg} or \code{NULL}.}    
    \item{dots}{ the list of  further arguments passed in \dots to be used in \code{\link{bootdist}} 
        in iterative calls to \code{\link{mledist}}, \code{\link{mmedist}},
        \code{\link{qmedist}}, \code{\link{mgedist}} or 
        \code{NULL} if no such arguments.}
    \item{convergence}{ an integer code for the convergence of 
    \code{\link{optim}}/\code{\link{constrOptim}} defined as below
    or defined by the user in the user-supplied optimization function. 
        \code{0} indicates successful convergence.
        \code{1} indicates that the iteration limit of \code{\link{optim}} has been reached.
        \code{10} indicates degeneracy of the Nealder-Mead simplex.
        \code{100} indicates that \code{\link{optim}} encountered an internal error.}    
    \item{discrete}{ the input argument or the automatic definition by the function to be passed
        to functions \code{\link{gofstat}}, \code{\link{plotdist}} 
        and \code{\link{cdfcomp}}. }
    
    \item{weights}{the vector of weigths used in the estimation process or \code{NULL}.}    

    Generic functions:
\describe{       
    \item{\code{print}}{
    The print of a \code{"fitdist"} object shows few traces about the fitting method and 
    the fitted distribution.
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
    \item{\code{logLik}}{ Extracts the estimated log-likelihood from the \code{"fitdist"} object.
    }
    \item{\code{vcov}}{ Extracts the estimated var-covariance matrix from the \code{"fitdist"} object
    (only available When \code{method = "mle"}).
    }
    \item{\code{coef}}{ Extracts the fitted coefficients from the \code{"fitdist"} object.
    }
}

}

\seealso{
    See \code{\link{fitdistrplus}} for an overview of the package.
    See \code{\link{mledist}}, \code{\link{mmedist}}, \code{\link{qmedist}}, 
    \code{\link{mgedist}} for details on parameter estimation.
    See \code{\link{gofstat}} for goodness-of-fit statistics.
    See \code{\link{plotdist}}, \code{\link{graphcomp}}, \code{\link{CIcdfplot}} for graphs
    (with or without uncertainty and/or multiple fits).
    See \code{\link{llplot}} for (log-)likelihood plots in the
    neighborhood of the fitted value.
    See \code{\link{bootdist}} for bootstrap procedures 
    and \code{\link{fitdistcens}} for censored-data fitting methods.
    See \code{\link{optim}} for base R optimization procedures.
    See \code{\link{quantile.fitdist}}, another generic function, which calculates 
    quantiles from the fitted distribution.
    See \code{\link{quantile}} for base R quantile computation.
}

\references{ 
Cullen AC and Frey HC (1999), \emph{Probabilistic techniques in exposure assessment}.
Plenum Press, USA, pp. 81-155.

Venables WN and Ripley BD (2002), \emph{Modern applied statistics with S}.
Springer, New York, pp. 435-446.

Vose D (2000), \emph{Risk analysis, a quantitative guide}.
John Wiley & Sons Ltd, Chischester, England, pp. 99-143.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}


\author{ 
Marie-Laure Delignette-Muller and Christophe Dutang.
}

\examples{

# (1) fit of a gamma distribution by maximum likelihood estimation
#

data(groundbeef)
serving <- groundbeef$serving
fitg <- fitdist(serving, "gamma")
summary(fitg)
plot(fitg)
plot(fitg, demp = TRUE)
plot(fitg, histo = FALSE, demp = TRUE)
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
gofstat(list(fitW, fitg, fitln), fitnames=c("Weibull", "gamma", "lognormal"))

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

cdfcomp(list(fitp,fitnb))
gofstat(list(fitp,fitnb))

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
res2 <- fitdist(mysample, "exp", start=mystart, custom.optim=myoptimize, 
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
    memp <- function(x, order) mean(x^order)

    #fit
    fP <- fitdist(x4, "pareto", method="mme", order=c(1, 2), memp="memp", 
      start=list(shape=10, scale=10), lower=1, upper=Inf)
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
cdfcomp(list(f1, f2, f3, f4, f5, f6, f7, f8), 
  xlogscale=TRUE, xlim=c(8, 250), verticals=TRUE)
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

# (15) Fit of a normal distribution on acute toxicity values of endosulfan in log10 for
# nonarthropod invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5 percent quantile value of 
# the fitted distribution (which is called the 5 percent hazardous concentration, HC5,
# in ecotoxicology) and estimation of other quantiles.
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

# (17) fit a non classical discrete distribution (the zero inflated Poisson distribution)
#
\dontrun{
require(gamlss.dist)
set.seed(1234)
x <- rZIP(n = 30, mu = 5, sigma = 0.2)
plotdist(x, discrete = TRUE)
fitzip <- fitdist(x, "ZIP", start =  list(mu = 4, sigma = 0.15), discrete = TRUE, 
  optim.method = "L-BFGS-B", lower = c(0, 0), upper = c(Inf, 1))
summary(fitzip)
plot(fitzip)
fitp <- fitdist(x, "pois")
cdfcomp(list(fitzip, fitp))
gofstat(list(fitzip, fitp))
}



# (18) examples with distributions in actuar (predefined starting values)
#
\dontrun{
require(actuar)
x <- c(2.3,0.1,2.7,2.2,0.4,2.6,0.2,1.,7.3,3.2,0.8,1.2,33.7,14.,
       21.4,7.7,1.,1.9,0.7,12.6,3.2,7.3,4.9,4000.,2.5,6.7,3.,63.,
       6.,1.6,10.1,1.2,1.5,1.2,30.,3.2,3.5,1.2,0.2,1.9,0.7,17.,
       2.8,4.8,1.3,3.7,0.2,1.8,2.6,5.9,2.6,6.3,1.4,0.8)
#log logistic
ft_llogis <- fitdist(x,'llogis')

x <- c(0.3837053, 0.8576858, 0.3552237, 0.6226119, 0.4783756, 0.3139799, 0.4051403, 
       0.4537631, 0.4711057, 0.5647414, 0.6479617, 0.7134207, 0.5259464, 0.5949068, 
       0.3509200, 0.3783077, 0.5226465, 1.0241043, 0.4384580, 1.3341520)
#inverse weibull
ft_iw <- fitdist(x,'invweibull')
}

}
\keyword{ distribution }
