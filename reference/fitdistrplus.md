# Overview of the fitdistrplus package

The idea of this package emerged in 2008 from a collaboration between
J.B. Denis, R. Pouillot and M.L. Delignette who at this time worked in
the area of quantitative risk assessment. The implementation of this
package was a part of a more general project named "Risk assessment with
R" gathering different packages and hosted in
[R-forge](https://r-forge.r-project.org/projects/riskassessment/).

The fitdistrplus package was first written by M.L. Delignette-Muller and
made available in
[CRAN](https://cran.r-project.org/package=fitdistrplus) on 2009 and
presented at the [2009 useR
conference](https://www.r-project.org/conferences/useR-2009/) in Rennes.
A few months after, C. Dutang joined the project by starting to
participate to the implementation of the fitdistrplus package. The
package has also been presented at the [2011 useR
conference](https://www.r-project.org/conferences/useR-2011/) at the
2eme rencontres R in 2013 (<https://r2013-lyon.sciencesconf.org/>), and
the [2019 useR
conference](https://www.r-project.org/conferences/useR-2019/). Since
2020, A. Siberchicot helps the development of fitdistrplus and maintains
the package.

Four vignettes are available within the package:

- a [general
  overview](https://lbbe-software.github.io/fitdistrplus/articles/fitdistrplus_vignette.html)
  of the package published in the Journal of Statistical Software
  ([doi:10.18637/jss.v064.i04](https://doi.org/10.18637/jss.v064.i04) ),

- a document answering the most [Frequently Asked
  Questions](https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html),

- a document presenting a [benchmark of optimization
  algorithms](https://lbbe-software.github.io/fitdistrplus/articles/Optimalgo.html)
  when finding parameters,

- a document about [starting
  values](https://lbbe-software.github.io/fitdistrplus/articles/fitdistrplus_vignette.html).

The fitdistrplus package is a general package that aims at helping the
fit of univariate parametric distributions to censored or non-censored
data. The two main functions are
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
for fit on non-censored data and
[`fitdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.md)
for fit on censored data.

The choice of candidate distributions to fit may be helped using
functions
[`descdist`](https://lbbe-software.github.io/fitdistrplus/reference/descdist.md)
and
[`plotdist`](https://lbbe-software.github.io/fitdistrplus/reference/plotdist.md)
for non-censored data and
[`plotdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/plotdistcens.md)
for censored data).

Using functions
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
and
[`fitdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.md),
different methods can be used to estimate the distribution parameters:

- maximum likelihood estimation by default
  ([`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)),

- moment matching estimation
  ([`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md)),

- quantile matching estimation
  ([`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md)),

- maximum goodness-of-fit estimation
  ([`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md)).

For classical distributions initial values are automatically calculated
if not provided by the user. Graphical functions
[`plotdist`](https://lbbe-software.github.io/fitdistrplus/reference/plotdist.md)
and
[`plotdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/plotdistcens.md)
can be used to help a manual calibration of initial values for
parameters of non-classical distributions. Function
[`prefit`](https://lbbe-software.github.io/fitdistrplus/reference/prefit.md)
is proposed to help the definition of good starting values in the
special case of constrained parameters. In the case where maximum
likelihood is chosen as the estimation method, function
[`llplot`](https://lbbe-software.github.io/fitdistrplus/reference/logLik-plot.md)
enables to visualize loglikelihood contours.

The goodness-of-fit of fitted distributions (a single fit or multiple
fits) can be explored using different graphical functions
([`cdfcomp`](https://lbbe-software.github.io/fitdistrplus/reference/graphcomp.md),
[`denscomp`](https://lbbe-software.github.io/fitdistrplus/reference/graphcomp.md),
[`qqcomp`](https://lbbe-software.github.io/fitdistrplus/reference/graphcomp.md)
and
[`ppcomp`](https://lbbe-software.github.io/fitdistrplus/reference/graphcomp.md)
for non-censored data and
[`cdfcompcens`](https://lbbe-software.github.io/fitdistrplus/reference/graphcompcens.md)
for censored data). Goodness-of-fit statistics are also provided for
non-censored data using function
[`gofstat`](https://lbbe-software.github.io/fitdistrplus/reference/gofstat.md).

Bootstrap is proposed to quantify the uncertainty on parameter estimates
(functions
[`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md)
and
[`bootdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/bootdistcens.md))
and also to quantify the uncertainty on CDF or quantiles estimated from
the fitted distribution
([`quantile`](https://lbbe-software.github.io/fitdistrplus/reference/quantile.md)
and
[`CIcdfplot`](https://lbbe-software.github.io/fitdistrplus/reference/CIcdfplot.md)).

## Author

Marie-Laure Delignette-Muller and Christophe Dutang.
