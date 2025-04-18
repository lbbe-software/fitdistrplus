\name{fitdistrplus-package}
\alias{fitdistrplus}
\alias{fitdistrplus-package}
\title{Overview of the \pkg{fitdistrplus} package}

\description{
The idea of this package emerged in 2008 from a collaboration between 
J.B. Denis, R. Pouillot and M.L. Delignette who at this time worked in the
area of quantitative risk assessment. The implementation of this package 
was a part of a more general project named "Risk assessment with R" 
gathering different packages and 
hosted in \href{https://r-forge.r-project.org/projects/riskassessment/}{R-forge}.

The \pkg{fitdistrplus} package was first written by M.L. Delignette-Muller and made 
available in \href{https://cran.r-project.org/package=fitdistrplus}{CRAN} on 2009 and presented at the \href{https://www.r-project.org/conferences/useR-2009/}{2009 useR conference} in Rennes.
A few months after, C. Dutang joined the project by
starting to participate to the implementation of the \pkg{fitdistrplus} package.
The package has also been presented at the \href{https://www.r-project.org/conferences/useR-2011/}{2011 useR conference}
at the 2eme rencontres R in 2013 (\url{https://r2013-lyon.sciencesconf.org/}),
and the \href{https://www.r-project.org/conferences/useR-2019/}{2019 useR conference}.
Since 2020, A. Siberchicot helps the development of \pkg{fitdistrplus} and maintains the package.

Four vignettes are available within the package:
\itemize{
\item a \href{https://lbbe-software.github.io/fitdistrplus/articles/fitdistrplus_vignette.html}{general overview} of the package published in the 
Journal of Statistical Software (\doi{10.18637/jss.v064.i04}),
\item a document answering the most \href{https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html}{Frequently Asked Questions},
\item a document presenting a \href{https://lbbe-software.github.io/fitdistrplus/articles/Optimalgo.html}{benchmark of optimization algorithms} when finding parameters,
\item a document about \href{https://lbbe-software.github.io/fitdistrplus/articles/fitdistrplus_vignette.html}{starting values}.
}

The \pkg{fitdistrplus} package is a general package that aims at helping the fit of univariate parametric 
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
distribution parameters: 
\itemize{
\item maximum likelihood estimation by default (\code{\link{mledist}}),
\item moment matching estimation (\code{\link{mmedist}}), 
\item quantile matching estimation (\code{\link{qmedist}}),
\item maximum goodness-of-fit estimation (\code{\link{mgedist}}). 
}
For classical distributions initial values are automatically calculated
if not provided by the user.
Graphical functions \code{\link{plotdist}} and \code{\link{plotdistcens}} 
can be used to help a manual calibration of initial values for parameters
of non-classical distributions. Function \code{\link{prefit}} is proposed
to help the definition of good starting values in the special case of
constrained parameters. In the case where maximum likelihood is chosen
as the estimation method, function \code{\link{llplot}} enables to
visualize loglikelihood contours.

The goodness-of-fit of fitted distributions (a single fit or multiple fits) can be explored
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

