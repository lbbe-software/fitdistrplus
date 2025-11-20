# Fit of univariate distributions to non-censored data

Fit of univariate distributions to non-censored data by maximum
likelihood (mle), moment matching (mme), quantile matching (qme) or
maximizing goodness-of-fit estimation (mge). The latter is also known as
minimizing distance estimation. Generic methods are `print`, `plot`,
`summary`, `quantile`, `logLik`, `AIC`, `BIC`, `vcov` and `coef`.

## Usage

``` r
fitdist(data, distr, method = c("mle", "mme", "qme", "mge", "mse"), 
    start=NULL, fix.arg=NULL, discrete, keepdata = TRUE, keepdata.nb=100, 
    calcvcov=TRUE, ...)
    
# S3 method for class 'fitdist'
print(x, ...)

# S3 method for class 'fitdist'
plot(x, breaks="default", ...)

# S3 method for class 'fitdist'
summary(object, ...)

# S3 method for class 'fitdist'
logLik(object, ...)

# S3 method for class 'fitdist'
AIC(object, ..., k = 2)

# S3 method for class 'fitdist'
BIC(object, ...)

# S3 method for class 'fitdist'
vcov(object, ...)

# S3 method for class 'fitdist'
coef(object, ...)
```

## Arguments

- data:

  A numeric vector.

- distr:

  A character string `"name"` naming a distribution for which the
  corresponding density function `dname`, the corresponding distribution
  function `pname` and the corresponding quantile function `qname` must
  be defined, or directly the density function.

- method:

  A character string coding for the fitting method: `"mle"` for 'maximum
  likelihood estimation', `"mme"` for 'moment matching estimation',
  `"qme"` for 'quantile matching estimation', `"mge"` for 'maximum
  goodness-of-fit estimation' and `"mse"` for 'maximum spacing
  estimation'.

- start:

  A named list giving the initial values of parameters of the named
  distribution or a function of data computing initial values and
  returning a named list. This argument may be omitted (default) for
  some distributions for which reasonable starting values are computed
  (see the 'details' section of
  [`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)).
  It may not be into account for closed-form formulas.

- fix.arg:

  An optional named list giving the values of fixed parameters of the
  named distribution or a function of data computing (fixed) parameter
  values and returning a named list. Parameters with fixed value are
  thus NOT estimated by this maximum likelihood procedure. The use of
  this argument is not possible if `method="mme"` and a closed-form
  formula is used.

- keepdata:

  a logical. If `TRUE`, dataset is returned, otherwise only a sample
  subset is returned.

- keepdata.nb:

  When `keepdata=FALSE`, the length (\>1) of the subset returned.

- calcvcov:

  A logical indicating if (asymptotic) covariance matrix is required.

- discrete:

  If TRUE, the distribution is considered as discrete. If `discrete` is
  missing, `discrete` is automaticaly set to `TRUE` when `distr` belongs
  to `"binom"`, `"nbinom"`, `"geom"`, `"hyper"` or `"pois"` and to
  `FALSE` in the other cases. It is thus recommended to enter this
  argument when using another discrete distribution. This argument will
  not directly affect the results of the fit but will be passed to
  functions
  [`gofstat`](https://lbbe-software.github.io/fitdistrplus/reference/gofstat.md),
  [`plotdist`](https://lbbe-software.github.io/fitdistrplus/reference/plotdist.md)
  and
  [`cdfcomp`](https://lbbe-software.github.io/fitdistrplus/reference/graphcomp.md).

- x:

  An object of class `"fitdist"`.

- object:

  An object of class `"fitdist"`.

- breaks:

  If `"default"` the histogram is plotted with the function `hist` with
  its default breaks definition. Else `breaks` is passed to the function
  `hist`. This argument is not taken into account with discrete
  distributions: `"binom"`, `"nbinom"`, `"geom"`, `"hyper"` and
  `"pois"`.

- k:

  penalty per parameter to be passed to the AIC generic function (2 by
  default).

- ...:

  Further arguments to be passed to generic functions, or to one of the
  functions `"mledist"`, `"mmedist"`, `"qmedist"` or `"mgedist"`
  depending of the chosen method. See
  [`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
  [`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md),
  [`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md),
  [`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md)
  for details on parameter estimation.

## Details

It is assumed that the `distr` argument specifies the distribution by
the probability density function, the cumulative distribution function
and the quantile function (d, p, q).

The four possible fitting methods are described below:

- When `method="mle"`:

  Maximum likelihood estimation consists in maximizing the
  log-likelihood. A numerical optimization is carried out in
  [`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)
  via `optim` to find the best values (see
  [`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)
  for details).

- When `method="mme"`:

  Moment matching estimation consists in equalizing theoretical and
  empirical moments. Estimated values of the distribution parameters are
  computed by a closed-form formula for the following distributions :
  `"norm"`, `"lnorm"`, `"pois"`, `"exp"`, `"gamma"`, `"nbinom"`,
  `"geom"`, `"beta"`, `"unif"` and `"logis"`. Otherwise the theoretical
  and the empirical moments are matched numerically, by minimization of
  the sum of squared differences between observed and theoretical
  moments. In this last case, further arguments are needed in the call
  to `fitdist`: `order` and `memp` (see
  [`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md)
  for details).

  Since Version 1.2-0, `mmedist` automatically computes the asymptotic
  covariance matrix, hence the theoretical moments `mdist` should be
  defined up to an order which equals to twice the maximal order given
  `order`.

- When `method = "qme"`:

  Quantile matching estimation consists in equalizing theoretical and
  empirical quantile. A numerical optimization is carried out in
  [`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md)
  via `optim` to minimize of the sum of squared differences between
  observed and theoretical quantiles. The use of this method requires an
  additional argument `probs`, defined as the numeric vector of the
  probabilities for which the quantile(s) is(are) to be matched (see
  [`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md)
  for details).

- When `method = "mge"`:

  Maximum goodness-of-fit estimation consists in maximizing a
  goodness-of-fit statistics. A numerical optimization is carried out in
  [`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md)
  via `optim` to minimize the goodness-of-fit distance. The use of this
  method requires an additional argument `gof` coding for the
  goodness-of-fit distance chosen. One can use the classical Cramer-von
  Mises distance (`"CvM"`), the classical Kolmogorov-Smirnov distance
  (`"KS"`), the classical Anderson-Darling distance (`"AD"`) which gives
  more weight to the tails of the distribution, or one of the variants
  of this last distance proposed by Luceno (2006) (see
  [`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md)
  for more details). This method is not suitable for discrete
  distributions.

- When `method = "mse"`:

  Maximum goodness-of-fit estimation consists in maximizing the average
  log spacing. A numerical optimization is carried out in
  [`msedist`](https://lbbe-software.github.io/fitdistrplus/reference/msedist.md)
  via `optim`.

By default, direct optimization of the log-likelihood (or other criteria
depending of the chosen method) is performed using
[`optim`](https://rdrr.io/r/stats/optim.html), with the "Nelder-Mead"
method for distributions characterized by more than one parameter and
the "BFGS" method for distributions characterized by only one parameter.
The optimization algorithm used in
[`optim`](https://rdrr.io/r/stats/optim.html) can be chosen or another
optimization function can be specified using ... argument (see
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)
for details). `start` may be omitted (i.e. `NULL`) for some classic
distributions (see the 'details' section of
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)).
Note that when errors are raised by `optim`, it's a good idea to start
by adding traces during the optimization process by adding
`control=list(trace=1, REPORT=1)` in ... argument.

Once the parameter(s) is(are) estimated, `fitdist` computes the
log-likelihood for every estimation method and for maximum likelihood
estimation the standard errors of the estimates calculated from the
Hessian at the solution found by `optim` or by the user-supplied
function passed to mledist.

By default (`keepdata = TRUE`), the object returned by `fitdist`
contains the data vector given in input. When dealing with large
datasets, we can remove the original dataset from the output by setting
`keepdata = FALSE`. In such a case, only `keepdata.nb` points (at most)
are kept by random subsampling `keepdata.nb`-2 points from the dataset
and adding the minimum and the maximum. If combined with
[`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md),
and use with non-parametric bootstrap be aware that bootstrap is
performed on the subset randomly selected in `fitdist`. Currently, the
graphical comparisons of multiple fits is not available in this
framework.

Weighted version of the estimation process is available for
`method = "mle", "mme", "qme"` by using `weights=...`. See the
corresponding man page for details. Weighted maximum GOF estimation
(when `method = "mge"`) is not allowed. It is not yet possible to take
into account weighths in functions `plotdist`, `plot.fitdist`,
`cdfcomp`, `denscomp`, `ppcomp`, `qqcomp`, `gofstat` and `descdist`
(developments planned in the future).

Once the parameter(s) is(are) estimated,
[`gofstat`](https://lbbe-software.github.io/fitdistrplus/reference/gofstat.md)
allows to compute goodness-of-fit statistics.

NB: if your data values are particularly small or large, a scaling may
be needed before the optimization process. See example (14) in this man
page and examples (14,15) in the test file of the package. Please also
take a look at the `Rmpfr` package available on CRAN for numerical
accuracy issues.

## Value

`fitdist` returns an object of class `"fitdist"`, a list with the
following components:

- estimate:

  the parameter estimates.

- method:

  the character string coding for the fitting method : `"mle"` for
  'maximum likelihood estimation', `"mme"` for 'matching moment
  estimation', `"qme"` for 'matching quantile estimation' `"mge"` for
  'maximum goodness-of-fit estimation' and `"mse"` for 'maximum spacing
  estimation'.

- sd:

  the estimated standard errors, `NA` if numerically not computable or
  `NULL` if not available.

- cor:

  the estimated correlation matrix, `NA` if numerically not computable
  or `NULL` if not available.

- vcov:

  the estimated variance-covariance matrix, `NULL` if not available for
  the estimation method considered.

- loglik:

  the log-likelihood.

- aic:

  the Akaike information criterion.

- bic:

  the the so-called BIC or SBC (Schwarz Bayesian criterion).

- n:

  the length of the data set.

- data:

  the data set.

- distname:

  the name of the distribution.

- fix.arg:

  the named list giving the values of parameters of the named
  distribution that must be kept fixed rather than estimated by maximum
  likelihood or `NULL` if there are no such parameters.

- fix.arg.fun:

  the function used to set the value of `fix.arg` or `NULL`.

- dots:

  the list of further arguments passed in ... to be used in
  [`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md)
  in iterative calls to
  [`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
  [`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md),
  [`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md),
  [`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md)
  or `NULL` if no such arguments.

- convergence:

  an integer code for the convergence of
  [`optim`](https://rdrr.io/r/stats/optim.html)/[`constrOptim`](https://rdrr.io/r/stats/constrOptim.html)
  defined as below or defined by the user in the user-supplied
  optimization function. `0` indicates successful convergence. `1`
  indicates that the iteration limit of
  [`optim`](https://rdrr.io/r/stats/optim.html) has been reached. `10`
  indicates degeneracy of the Nealder-Mead simplex. `100` indicates that
  [`optim`](https://rdrr.io/r/stats/optim.html) encountered an internal
  error.

- discrete:

  the input argument or the automatic definition by the function to be
  passed to functions
  [`gofstat`](https://lbbe-software.github.io/fitdistrplus/reference/gofstat.md),
  [`plotdist`](https://lbbe-software.github.io/fitdistrplus/reference/plotdist.md)
  and
  [`cdfcomp`](https://lbbe-software.github.io/fitdistrplus/reference/graphcomp.md).

- weights:

  the vector of weigths used in the estimation process or `NULL`.

Generic functions:

- `print`:

  The print of a `"fitdist"` object shows few traces about the fitting
  method and the fitted distribution.

- `summary`:

  The summary provides the parameter estimates of the fitted
  distribution, the log-likelihood, AIC and BIC statistics and when the
  maximum likelihood is used, the standard errors of the parameter
  estimates and the correlation matrix between parameter estimates.

- `plot`:

  The plot of an object of class "fitdist" returned by `fitdist` uses
  the function
  [`plotdist`](https://lbbe-software.github.io/fitdistrplus/reference/plotdist.md).
  An object of class "fitdist" or a list of objects of class "fitdist"
  corresponding to various fits using the same data set may also be
  plotted using a cdf plot (function
  [`cdfcomp`](https://lbbe-software.github.io/fitdistrplus/reference/graphcomp.md)),
  a density plot(function
  [`denscomp`](https://lbbe-software.github.io/fitdistrplus/reference/graphcomp.md)),
  a density Q-Q plot (function
  [`qqcomp`](https://lbbe-software.github.io/fitdistrplus/reference/graphcomp.md)),
  or a P-P plot (function
  [`ppcomp`](https://lbbe-software.github.io/fitdistrplus/reference/graphcomp.md)).

- `logLik`:

  Extracts the estimated log-likelihood from the `"fitdist"` object.

- `AIC`:

  Extracts the AIC from the `"fitdist"` object.

- `BIC`:

  Extracts the estimated BIC from the `"fitdist"` object.

- `vcov`:

  Extracts the estimated var-covariance matrix from the `"fitdist"`
  object (only available When `method = "mle"`).

- `coef`:

  Extracts the fitted coefficients from the `"fitdist"` object.

## See also

See
[`fitdistrplus`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistrplus.md)
for an overview of the package. See
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
[`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md),
[`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md),
[`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md),
[`msedist`](https://lbbe-software.github.io/fitdistrplus/reference/msedist.md)
for details on parameter estimation. See
[`gofstat`](https://lbbe-software.github.io/fitdistrplus/reference/gofstat.md)
for goodness-of-fit statistics. See
[`plotdist`](https://lbbe-software.github.io/fitdistrplus/reference/plotdist.md),
[`graphcomp`](https://lbbe-software.github.io/fitdistrplus/reference/graphcomp.md),
[`CIcdfplot`](https://lbbe-software.github.io/fitdistrplus/reference/CIcdfplot.md)
for graphs (with or without uncertainty and/or multiple fits). See
[`llplot`](https://lbbe-software.github.io/fitdistrplus/reference/logLik-plot.md)
for (log-)likelihood plots in the neighborhood of the fitted value. See
[`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md)
for bootstrap procedures and
[`fitdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.md)
for censored-data fitting methods. See
[`optim`](https://rdrr.io/r/stats/optim.html) for base R optimization
procedures. See
[`quantile.fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/quantile.md),
another generic function, which calculates quantiles from the fitted
distribution. See
[`quantile`](https://lbbe-software.github.io/fitdistrplus/reference/quantile.md)
for base R quantile computation.

Please visit the [Frequently Asked
Questions](https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html).

## References

I. Ibragimov and R. Has'minskii (1981), *Statistical Estimation -
Asymptotic Theory*, Springer-Verlag,
[doi:10.1007/978-1-4899-0027-2](https://doi.org/10.1007/978-1-4899-0027-2)

Cullen AC and Frey HC (1999), *Probabilistic techniques in exposure
assessment*. Plenum Press, USA, pp. 81-155.

Venables WN and Ripley BD (2002), *Modern applied statistics with S*.
Springer, New York, pp. 435-446,
[doi:10.1007/978-0-387-21706-2](https://doi.org/10.1007/978-0-387-21706-2)
.

Vose D (2000), *Risk analysis, a quantitative guide*. John Wiley & Sons
Ltd, Chischester, England, pp. 99-143.

Delignette-Muller ML and Dutang C (2015), *fitdistrplus: An R Package
for Fitting Distributions*. Journal of Statistical Software, 64(4),
1-34, [doi:10.18637/jss.v064.i04](https://doi.org/10.18637/jss.v064.i04)
.

## Author

Marie-Laure Delignette-Muller and Christophe Dutang.

## Examples

``` r
# (1) fit of a gamma distribution by maximum likelihood estimation
#

data(groundbeef)
serving <- groundbeef$serving
fitg <- fitdist(serving, "gamma")
summary(fitg)
#> Fitting of the distribution ' gamma ' by maximum likelihood 
#> Parameters : 
#>         estimate  Std. Error
#> shape 4.00955898 0.341451640
#> rate  0.05443907 0.004937239
#> Loglikelihood:  -1253.625   AIC:  2511.25   BIC:  2518.325 
#> Correlation matrix:
#>           shape      rate
#> shape 1.0000000 0.9384578
#> rate  0.9384578 1.0000000
#> 
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
#> Fitting of the distribution ' gamma ' by matching moments 
#> Parameters : 
#>         estimate  Std. Error
#> shape 4.22848617 0.417232914
#> rate  0.05741663 0.005930118
#> Loglikelihood:  -1253.825   AIC:  2511.65   BIC:  2518.724 
#> Correlation matrix:
#>           shape      rate
#> shape 1.0000000 0.9553622
#> rate  0.9553622 1.0000000
#> 

# (3) Comparison of various fits 
#

fitW <- fitdist(serving, "weibull")
fitg <- fitdist(serving, "gamma")
fitln <- fitdist(serving, "lnorm")
summary(fitW)
#> Fitting of the distribution ' weibull ' by maximum likelihood 
#> Parameters : 
#>        estimate Std. Error
#> shape  2.185885  0.1045755
#> scale 83.347679  2.5268631
#> Loglikelihood:  -1255.225   AIC:  2514.449   BIC:  2521.524 
#> Correlation matrix:
#>          shape    scale
#> shape 1.000000 0.321821
#> scale 0.321821 1.000000
#> 
summary(fitg)
#> Fitting of the distribution ' gamma ' by maximum likelihood 
#> Parameters : 
#>         estimate  Std. Error
#> shape 4.00955898 0.341451640
#> rate  0.05443907 0.004937239
#> Loglikelihood:  -1253.625   AIC:  2511.25   BIC:  2518.325 
#> Correlation matrix:
#>           shape      rate
#> shape 1.0000000 0.9384578
#> rate  0.9384578 1.0000000
#> 
summary(fitln)
#> Fitting of the distribution ' lnorm ' by maximum likelihood 
#> Parameters : 
#>          estimate Std. Error
#> meanlog 4.1693701 0.03366988
#> sdlog   0.5366095 0.02380783
#> Loglikelihood:  -1261.319   AIC:  2526.639   BIC:  2533.713 
#> Correlation matrix:
#>         meanlog sdlog
#> meanlog       1     0
#> sdlog         0     1
#> 
cdfcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))

denscomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))

qqcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))

ppcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))

gofstat(list(fitW, fitg, fitln), fitnames=c("Weibull", "gamma", "lognormal"))
#> Goodness-of-fit statistics
#>                                Weibull     gamma lognormal
#> Kolmogorov-Smirnov statistic 0.1396646 0.1281486 0.1493090
#> Cramer-von Mises statistic   0.6840994 0.6936274 0.8277358
#> Anderson-Darling statistic   3.5736460 3.5672625 4.5436542
#> 
#> Goodness-of-fit criteria
#>                                 Weibull    gamma lognormal
#> Akaike's Information Criterion 2514.449 2511.250  2526.639
#> Bayesian Information Criterion 2521.524 2518.325  2533.713

# (4) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view 
# dedicated to probability distributions
#

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

fitgumbel <- fitdist(serving, "gumbel", start=list(a=10, b=10))
#> Error in fitdist(serving, "gumbel", start = list(a = 10, b = 10)): The  dgumbel  function must be defined
summary(fitgumbel)
#> Error: object 'fitgumbel' not found
plot(fitgumbel)
#> Error: object 'fitgumbel' not found

# (5) fit discrete distributions (Poisson and negative binomial)
#

data(toxocara)
number <- toxocara$number
fitp <- fitdist(number,"pois")
summary(fitp)
#> Fitting of the distribution ' pois ' by maximum likelihood 
#> Parameters : 
#>        estimate Std. Error
#> lambda 8.679245  0.4046719
#> Loglikelihood:  -507.5334   AIC:  1017.067   BIC:  1019.037 
plot(fitp)

fitnb <- fitdist(number,"nbinom")
summary(fitnb)
#> Fitting of the distribution ' nbinom ' by maximum likelihood 
#> Parameters : 
#>       estimate Std. Error
#> size 0.3971457 0.08289027
#> mu   8.6802520 1.93501002
#> Loglikelihood:  -159.3441   AIC:  322.6882   BIC:  326.6288 
#> Correlation matrix:
#>              size           mu
#> size  1.000000000 -0.000103854
#> mu   -0.000103854  1.000000000
#> 
plot(fitnb)


cdfcomp(list(fitp,fitnb))

gofstat(list(fitp,fitnb))
#> Chi-squared statistic:  31256.96 7.48606 
#> Degree of freedom of the Chi-squared distribution:  5 4 
#> Chi-squared p-value:  0 0.1123255 
#>    the p-value may be wrong with some theoretical counts < 5  
#> Chi-squared table:
#>       obscounts theo 1-mle-pois theo 2-mle-nbinom
#> <= 0         14     0.009014207         15.295027
#> <= 1          8     0.078236515          5.808596
#> <= 3          6     1.321767253          6.845015
#> <= 4          6     2.131297825          2.407815
#> <= 9          6    29.827829425          7.835196
#> <= 21         6    19.626223437          8.271110
#> > 21          7     0.005631338          6.537242
#> 
#> Goodness-of-fit criteria
#>                                1-mle-pois 2-mle-nbinom
#> Akaike's Information Criterion   1017.067     322.6882
#> Bayesian Information Criterion   1019.037     326.6288

# (6) how to change the optimisation method?
#

data(groundbeef)
serving <- groundbeef$serving
fitdist(serving, "gamma", optim.method="Nelder-Mead")
#> Fitting of the distribution ' gamma ' by maximum likelihood 
#> Parameters:
#>         estimate  Std. Error
#> shape 4.00955898 0.341451640
#> rate  0.05443907 0.004937239
fitdist(serving, "gamma", optim.method="BFGS") 
#> Fitting of the distribution ' gamma ' by maximum likelihood 
#> Parameters:
#>         estimate  Std. Error
#> shape 4.21183650 0.359345675
#> rate  0.05719298 0.005180917
fitdist(serving, "gamma", optim.method="SANN")
#> Fitting of the distribution ' gamma ' by maximum likelihood 
#> Parameters:
#>         estimate  Std. Error
#> shape 3.94902092 0.336099391
#> rate  0.05370886 0.004872886

# (7) custom optimization function
#
# \donttest{
#create the sample
set.seed(1234)
mysample <- rexp(100, 5)
mystart <- list(rate=8)

res1 <- fitdist(mysample, dexp, start= mystart, optim.method="Nelder-Mead")

#show the result
summary(res1)
#> Fitting of the distribution ' exp ' by maximum likelihood 
#> Parameters : 
#>      estimate Std. Error
#> rate 5.120312  0.5120312
#> Loglikelihood:  63.32596   AIC:  -124.6519   BIC:  -122.0467 

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
#> Fitting of the distribution ' exp ' by maximum likelihood 
#> Parameters : 
#>      estimate
#> rate 5.120531
#> Loglikelihood:  63.32596   AIC:  -124.6519   BIC:  -122.0467 
# }


# (8) custom optimization function - another example with the genetic algorithm
#
# \donttest{
    #set a sample
    fit1 <- fitdist(serving, "gamma")
    summary(fit1)
#> Fitting of the distribution ' gamma ' by maximum likelihood 
#> Parameters : 
#>         estimate  Std. Error
#> shape 4.00955898 0.341451640
#> rate  0.05443907 0.004937239
#> Loglikelihood:  -1253.625   AIC:  2511.25   BIC:  2518.325 
#> Correlation matrix:
#>           shape      rate
#> shape 1.0000000 0.9384578
#> rate  0.9384578 1.0000000
#> 

    #wrap genoud function rgenoud package
    mygenoud <- function(fn, par, ...) 
    {
        require("rgenoud")
        res <- genoud(fn, starting.values=par, ...)        
        standardres <- c(res, convergence=0)
            
        return(standardres)
    }

    #call fitdist with a 'custom' optimization function
    fit2 <- fitdist(serving, "gamma", custom.optim=mygenoud, nvars=2,    
        Domains=cbind(c(0, 0), c(10, 10)), boundary.enforcement=1, 
        print.level=1, hessian=TRUE)
#> Loading required package: rgenoud
#> ##  rgenoud (Version 5.9-0.11, Build Date: 2024-10-03)
#> ##  See http://sekhon.berkeley.edu/rgenoud for additional documentation.
#> ##  Please cite software as:
#> ##   Walter Mebane, Jr. and Jasjeet S. Sekhon. 2011.
#> ##   ``Genetic Optimization Using Derivatives: The rgenoud package for R.''
#> ##   Journal of Statistical Software, 42(11): 1-26. 
#> ##
#> 
#> 
#> Thu Nov 20 20:46:21 2025
#> Domains:
#>  0.000000e+00   <=  X1   <=    1.000000e+01 
#>  0.000000e+00   <=  X2   <=    1.000000e+01 
#> 
#> Data Type: Floating Point
#> Operators (code number, name, population) 
#>  (1) Cloning...........................  122
#>  (2) Uniform Mutation..................  125
#>  (3) Boundary Mutation.................  125
#>  (4) Non-Uniform Mutation..............  125
#>  (5) Polytope Crossover................  125
#>  (6) Simple Crossover..................  126
#>  (7) Whole Non-Uniform Mutation........  125
#>  (8) Heuristic Crossover...............  126
#>  (9) Local-Minimum Crossover...........  0
#> 
#> HARD Maximum Number of Generations: 100
#> Maximum Nonchanging Generations: 10
#> Population size       : 1000
#> Convergence Tolerance: 1.000000e-03
#> 
#> Using the BFGS Derivative Based Optimizer on the Best Individual Each Generation.
#> Checking Gradients before Stopping.
#> Not Using Out of Bounds Individuals But Allowing Trespassing.
#> 
#> Minimization Problem.
#> 
#> 
#> Generation#      Solution Value
#> 
#>       0  4.936206e+00
#> 
#> 'wait.generations' limit reached.
#> No significant improvement in 10 generations.
#> 
#> Solution Fitness Value: 4.935532e+00
#> 
#> Parameters at the Solution (parameter, gradient):
#> 
#>  X[ 1] : 4.008339e+00    G[ 1] : 6.842735e-09
#>  X[ 2] : 5.442736e-02    G[ 2] : -1.863593e-07
#> 
#> Solution Found Generation 1
#> Number of Generations Run 11
#> 
#> Thu Nov 20 20:46:22 2025
#> Total run time : 0 hours 0 minutes and 1 seconds

    summary(fit2)
#> Fitting of the distribution ' gamma ' by maximum likelihood 
#> Parameters : 
#>         estimate  Std. Error
#> shape 4.00833917 0.341343852
#> rate  0.05442736 0.004936216
#> Loglikelihood:  -1253.625   AIC:  2511.25   BIC:  2518.325 
#> Correlation matrix:
#>           shape      rate
#> shape 1.0000000 0.9384395
#> rate  0.9384395 1.0000000
#> 
# }

# (9) estimation of the standard deviation of a gamma distribution 
# by maximum likelihood with the shape fixed at 4 using the argument fix.arg
#

data(groundbeef)
serving <- groundbeef$serving
f1c  <- fitdist(serving,"gamma",start=list(rate=0.1),fix.arg=list(shape=4))
summary(f1c)
#> Fitting of the distribution ' gamma ' by maximum likelihood 
#> Parameters : 
#>        estimate  Std. Error
#> rate 0.05431619 0.001703472
#> Fixed parameters:
#>       value
#> shape     4
#> Loglikelihood:  -1253.625   AIC:  2509.251   BIC:  2512.788 
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
#> Fitting of the distribution ' weibull ' by maximum likelihood 
#> Parameters : 
#>        estimate Std. Error
#> shape  2.185885  0.1045755
#> scale 83.347679  2.5268631
#> Loglikelihood:  -1255.225   AIC:  2514.449   BIC:  2521.524 
#> Correlation matrix:
#>          shape    scale
#> shape 1.000000 0.321821
#> scale 0.321821 1.000000
#> 
plot(fWmle)

gofstat(fWmle)
#> Goodness-of-fit statistics
#>                              1-mle-weibull
#> Kolmogorov-Smirnov statistic     0.1396646
#> Cramer-von Mises statistic       0.6840994
#> Anderson-Darling statistic       3.5736460
#> 
#> Goodness-of-fit criteria
#>                                1-mle-weibull
#> Akaike's Information Criterion      2514.449
#> Bayesian Information Criterion      2521.524

fWqme <- fitdist(serving, "weibull", method="qme", probs=c(0.25, 0.75))
summary(fWqme)
#> Fitting of the distribution ' weibull ' by matching quantiles 
#> Parameters : 
#>        estimate
#> shape  2.268699
#> scale 86.590853
#> Loglikelihood:  -1256.129   AIC:  2516.258   BIC:  2523.332 
plot(fWqme)

gofstat(fWqme)
#> Goodness-of-fit statistics
#>                              1-qme-weibull
#> Kolmogorov-Smirnov statistic     0.1692858
#> Cramer-von Mises statistic       0.9664709
#> Anderson-Darling statistic       4.8479858
#> 
#> Goodness-of-fit criteria
#>                                1-qme-weibull
#> Akaike's Information Criterion      2516.258
#> Bayesian Information Criterion      2523.332


# (11) Fit of a Pareto distribution by numerical moment matching estimation
#

# \donttest{
    require("actuar")
#> Loading required package: actuar
#> 
#> Attaching package: ‘actuar’
#> The following objects are masked from ‘package:stats’:
#> 
#>     sd, var
#> The following object is masked from ‘package:grDevices’:
#> 
#>     cm
    #simulate a sample
    x4 <- rpareto(1000, 6, 2)

    #empirical raw moment
    memp <- function(x, order) mean(x^order)

    #fit
    fP <- fitdist(x4, "pareto", method="mme", order=c(1, 2), memp="memp", 
      start=list(shape=10, scale=10), lower=1, upper=Inf)
#> Error in mmedist(data, distname, start = arg_startfix$start.arg, fix.arg = arg_startfix$fix.arg,     checkstartfix = TRUE, calcvcov = calcvcov, ...): the empirical moment must be defined as a function
    summary(fP)
#> Error: object 'fP' not found
    plot(fP)
#> Error: object 'fP' not found

# }

# (12) Fit of a Weibull distribution to serving size data by maximum 
# goodness-of-fit estimation using all the distances available
# 
# \donttest{
data(groundbeef)
serving <- groundbeef$serving
(f1 <- fitdist(serving, "weibull", method="mge", gof="CvM"))
#> Fitting of the distribution ' weibull ' by maximum goodness-of-fit 
#> Parameters:
#>        estimate
#> shape  2.093204
#> scale 82.660014
(f2 <- fitdist(serving, "weibull", method="mge", gof="KS"))
#> Fitting of the distribution ' weibull ' by maximum goodness-of-fit 
#> Parameters:
#>        estimate
#> shape  2.065634
#> scale 81.450487
(f3 <- fitdist(serving, "weibull", method="mge", gof="AD"))
#> Fitting of the distribution ' weibull ' by maximum goodness-of-fit 
#> Parameters:
#>        estimate
#> shape  2.125425
#> scale 82.890502
(f4 <- fitdist(serving, "weibull", method="mge", gof="ADR"))
#> Fitting of the distribution ' weibull ' by maximum goodness-of-fit 
#> Parameters:
#>        estimate
#> shape  2.072035
#> scale 82.762593
(f5 <- fitdist(serving, "weibull", method="mge", gof="ADL"))
#> Fitting of the distribution ' weibull ' by maximum goodness-of-fit 
#> Parameters:
#>        estimate
#> shape  2.197498
#> scale 82.016005
(f6 <- fitdist(serving, "weibull", method="mge", gof="AD2R"))
#> Fitting of the distribution ' weibull ' by maximum goodness-of-fit 
#> Parameters:
#>       estimate
#> shape  1.90328
#> scale 81.33464
(f7 <- fitdist(serving, "weibull", method="mge", gof="AD2L"))
#> Fitting of the distribution ' weibull ' by maximum goodness-of-fit 
#> Parameters:
#>        estimate
#> shape  2.483836
#> scale 78.252113
(f8 <- fitdist(serving, "weibull", method="mge", gof="AD2"))
#> Fitting of the distribution ' weibull ' by maximum goodness-of-fit 
#> Parameters:
#>        estimate
#> shape  2.081168
#> scale 85.281194
cdfcomp(list(f1, f2, f3, f4, f5, f6, f7, f8))

cdfcomp(list(f1, f2, f3, f4, f5, f6, f7, f8), 
  xlogscale=TRUE, xlim=c(8, 250), verticals=TRUE)

denscomp(list(f1, f2, f3, f4, f5, f6, f7, f8))

# }

# (13) Fit of a uniform distribution using maximum likelihood 
# (a closed formula is used in this special case where the loglikelihood is not defined),
# or maximum goodness-of-fit with Cramer-von Mises or Kolmogorov-Smirnov distance
# 

set.seed(1234)
u <- runif(50, min=5, max=10)

fumle <- fitdist(u, "unif", method="mle")
summary(fumle)
#> Fitting of the distribution ' unif ' by maximum likelihood 
#> Parameters : 
#>     estimate
#> min 5.047479
#> max 9.960752
#> Loglikelihood:  -79.59702   AIC:  163.194   BIC:  167.0181 
plot(fumle)

gofstat(fumle)
#> Goodness-of-fit statistics
#>                              1-mle-unif
#> Kolmogorov-Smirnov statistic  0.1340723
#> Cramer-von Mises statistic    0.1566892
#> Anderson-Darling statistic          Inf
#> 
#> Goodness-of-fit criteria
#>                                1-mle-unif
#> Akaike's Information Criterion   163.1940
#> Bayesian Information Criterion   167.0181

fuCvM <- fitdist(u, "unif", method="mge", gof="CvM")
summary(fuCvM)
#> Fitting of the distribution ' unif ' by maximum goodness-of-fit 
#> Parameters : 
#>     estimate
#> min 5.110497
#> max 9.552878
#> Loglikelihood:  -Inf   AIC:  Inf   BIC:  Inf 
plot(fuCvM)

gofstat(fuCvM)
#> Goodness-of-fit statistics
#>                              1-mge-unif
#> Kolmogorov-Smirnov statistic 0.11370966
#> Cramer-von Mises statistic   0.07791651
#> Anderson-Darling statistic          Inf
#> 
#> Goodness-of-fit criteria
#>                                1-mge-unif
#> Akaike's Information Criterion        Inf
#> Bayesian Information Criterion        Inf

fuKS <- fitdist(u, "unif", method="mge", gof="KS")
summary(fuKS)
#> Fitting of the distribution ' unif ' by maximum goodness-of-fit 
#> Parameters : 
#>     estimate
#> min 5.092357
#> max 9.323818
#> Loglikelihood:  -Inf   AIC:  Inf   BIC:  Inf 
plot(fuKS)

gofstat(fuKS)
#> Goodness-of-fit statistics
#>                              1-mge-unif
#> Kolmogorov-Smirnov statistic 0.09216159
#> Cramer-von Mises statistic   0.12241830
#> Anderson-Darling statistic          Inf
#> 
#> Goodness-of-fit criteria
#>                                1-mge-unif
#> Akaike's Information Criterion        Inf
#> Bayesian Information Criterion        Inf

# (14) scaling problem
# the simulated dataset (below) has particularly small values, hence without scaling (10^0),
# the optimization raises an error. The for loop shows how scaling by 10^i
# for i=1,...,6 makes the fitting procedure work correctly.

set.seed(1234)
x2 <- rnorm(100, 1e-4, 2e-4)

for(i in 0:6)
        cat(i, try(fitdist(x2*10^i, "cauchy", method="mle")$estimate, silent=TRUE), "\n")
#> <simpleError in optim(par = vstart, fn = fnobj, fix.arg = fix.arg, obs = data,     gr = gradient, ddistnam = ddistname, hessian = TRUE, method = meth,     lower = lower, upper = upper, ...): non-finite finite-difference value [2]>
#> 0 Error in fitdist(x2 * 10^i, "cauchy", method = "mle") : 
#>   the function mle failed to estimate the parameters, 
#>                 with the error code 100
#> 
#>  
#> <simpleError in optim(par = vstart, fn = fnobj, fix.arg = fix.arg, obs = data,     gr = gradient, ddistnam = ddistname, hessian = TRUE, method = meth,     lower = lower, upper = upper, ...): non-finite finite-difference value [2]>
#> 1 Error in fitdist(x2 * 10^i, "cauchy", method = "mle") : 
#>   the function mle failed to estimate the parameters, 
#>                 with the error code 100
#> 
#>  
#> 2 0.001870693 0.01100646 
#> 3 0.01871473 0.1100713 
#> 4 0.1870693 1.100646 
#> 5 1.876032 11.0131 
#> 6 18.76032 110.131 

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
#> Estimated quantiles for each specified probability (non-censored data)
#>            p=0.05
#> estimate 1.744227
quantile(fln, probs = c(0.05, 0.1, 0.2))
#> Estimated quantiles for each specified probability (non-censored data)
#>            p=0.05    p=0.1  p=0.2
#> estimate 1.744227 2.080093 2.4868


# (16) Fit of a triangular distribution using Cramer-von Mises or
# Kolmogorov-Smirnov distance
# 

# \donttest{
set.seed(1234)
require("mc2d")
#> Loading required package: mc2d
#> Loading required package: mvtnorm
#> 
#> Attaching package: ‘mc2d’
#> The following objects are masked from ‘package:base’:
#> 
#>     pmax, pmin
t <- rtriang(100, min=5, mode=6, max=10)
fCvM <- fitdist(t, "triang", method="mge", start = list(min=4, mode=6,max=9), gof="CvM")
#> Warning: Some parameter names have no starting/fixed value but have a default value: mean.
fKS <- fitdist(t, "triang", method="mge", start = list(min=4, mode=6,max=9), gof="KS")
#> Warning: Some parameter names have no starting/fixed value but have a default value: mean.
cdfcomp(list(fCvM,fKS))

# }

# (17) fit a non classical discrete distribution (the zero inflated Poisson distribution)
#
# \donttest{
require("gamlss.dist")
#> Loading required package: gamlss.dist
set.seed(1234)
x <- rZIP(n = 30, mu = 5, sigma = 0.2)
plotdist(x, discrete = TRUE)

fitzip <- fitdist(x, "ZIP", start =  list(mu = 4, sigma = 0.15), discrete = TRUE, 
  optim.method = "L-BFGS-B", lower = c(0, 0), upper = c(Inf, 1))
#> Warning: The dZIP function should return a zero-length vector when input has length zero
#> Warning: The pZIP function should return a zero-length vector when input has length zero
summary(fitzip)
#> Fitting of the distribution ' ZIP ' by maximum likelihood 
#> Parameters : 
#>        estimate Std. Error
#> mu    4.3166098 0.43412153
#> sigma 0.1891794 0.07416889
#> Loglikelihood:  -67.13886   AIC:  138.2777   BIC:  141.0801 
#> Correlation matrix:
#>               mu      sigma
#> mu    1.00000000 0.06418931
#> sigma 0.06418931 1.00000000
#> 
plot(fitzip)

fitp <- fitdist(x, "pois")
cdfcomp(list(fitzip, fitp))

gofstat(list(fitzip, fitp))
#> Chi-squared statistic:  3.579708 35.91516 
#> Degree of freedom of the Chi-squared distribution:  3 4 
#> Chi-squared p-value:  0.3105704 3.012341e-07 
#>    the p-value may be wrong with some theoretical counts < 5  
#> Chi-squared table:
#>      obscounts theo 1-mle-ZIP theo 2-mle-pois
#> <= 0         6       5.999996       0.9059215
#> <= 2         7       4.425507       8.7194943
#> <= 4         5       9.047522      12.1379326
#> <= 5         5       4.054142       3.9650580
#> <= 7         5       4.715294       3.4694258
#> > 7          2       1.757539       0.8021677
#> 
#> Goodness-of-fit criteria
#>                                1-mle-ZIP 2-mle-pois
#> Akaike's Information Criterion  138.2777   153.7397
#> Bayesian Information Criterion  141.0801   155.1409
# }



# (18) examples with distributions in actuar (predefined starting values)
#
# \donttest{
require("actuar")
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
# }
```
