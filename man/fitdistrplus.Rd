\name{fitdistrplus-package}
\alias{fitdistrplus}
\title{General presentation of the package}

\description{
The idea of this package emerged in 2008 from a collaboration between 
JB Denis, R Pouillot and ML Delignette who at this time worked in the
area of quantitative risk assessment. The implementation of this package 
was a part of a more general project named "Risk assessment with R" 
gathering different packages and 
hosted in Rforge (\url{http://r-forge.r-project.org/projects/riskassessment/}).
The fitdistrplus package was first written by ML Delignette-Muller and made 
available in CRAN on 2009.
A few months after C Dutang joined the project by
starting to participate to the implementation of the fitdistrplus package.

Two vignettes are available within the package. The first one corresponds to 
a general presentation of the package published in the Journal of Statistical Software
and the second is a document answering the most Frequently Asked Questions (a FAQ vignette).

Fitdistrplus is a general package that aims at helping the fit of univariate parametric 
distributions to censored or
non-censored data. The two main functions are 
\code{\link{fitdist}} for fit on non-censored data and 
\code{\link{fitdistcens}} for fit on censored data. 

The choice of candidate 
distributions to fit may be helped using functions \code{\link{descdist}} and
\code{\link{plotdist}} for non-censored data and \code{\link{plotdistcens}}
for censored data). 

Using functions \code{\link{fitdist}} and 
\code{\link{fitdistcens}}, different methods can be used to estimate the 
distribution parameters: maximum likelihood by default (\code{\link{mledist}}),
moment matching (\code{\link{mmedist}}), quantile matching (\code{\link{qmedist}})
or maximum goodness-of-fit estimation (\code{\link{mgedist}}). For classical
distributions initial values are automatically calculated if not provided by the user.
Graphical functions \code{\link{plotdist}} and \code{\link{plotdistcens}} 
can be used to help a manual calibration of initial values for parameters
of non-classical distributions. Function \code{\link{prefit}} is proposed
to help the definition of good starting values in the special case of
constrained parameters. In the case where maximum likelihood is chosen
as the estimation method, function \code{\link{llplot}} enables to
visualize loglikelihood contours.

The goodness-of-fit of fitted distributions can be explored
using different graphical functions (\code{\link{cdfcomp}}, \code{\link{denscomp}},
\code{\link{qqcomp}} and \code{\link{ppcomp}} for non-censored data and
\code{\link{cdfcompcens}} for censored data). Goodness-of-fit statistics are also 
provided for non-censored data using function \code{\link{gofstat}}.

Bootstrap is proposed to quantify the uncertainty on parameter estimates
(functions \code{\link{bootdist}} and \code{\link{bootdistcens}})
and also to quantify the uncertainty on CDF or quantiles estimated
from the fitted distribution (\code{\link{quantile}} and \code{\link{CIcdfplot}}).
}

\author{ 
Marie-Laure Delignette-Muller and Christophe Dutang.
}

