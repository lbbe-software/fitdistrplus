# Bootstrap simulation of uncertainty for non-censored data

Uses parametric or nonparametric bootstrap resampling in order to
simulate uncertainty in the parameters of the distribution fitted to
non-censored data.

## Usage

``` r
bootdist(f, bootmethod = "param", niter = 1001, silent = TRUE, 
      parallel = c("no", "snow", "multicore"), ncpus)
# S3 method for class 'bootdist'
print(x, ...)
# S3 method for class 'bootdist'
plot(x, main = "Bootstrapped values of parameters", enhance = FALSE, 
    trueval = NULL, rampcol = NULL, nbgrid = 100, nbcol = 100, ...)
# S3 method for class 'bootdist'
summary(object, ...)
# S3 method for class 'bootdist'
density(..., bw = nrd0, adjust = 1, kernel = "gaussian")
# S3 method for class 'density.bootdist'
plot(x, mar=c(4,4,2,1), lty=NULL, col=NULL, lwd=NULL, ...)
# S3 method for class 'density.bootdist'
print(x, ...)
```

## Arguments

- f:

  An object of class `"fitdist"`, output of the
  [`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
  function.

- bootmethod:

  A character string coding for the type of resampling : `"param"` for a
  parametric resampling and `"nonparam"` for a nonparametric resampling
  of data.

- niter:

  The number of samples drawn by bootstrap.

- silent:

  A logical to remove or show warnings and errors when bootstraping.

- parallel:

  The type of parallel operation to be used, `"snow"` or `"multicore"`
  (the second one not being available on Windows), or `"no"` if no
  parallel operation.

- ncpus:

  Number of processes to be used in parallel operation : typically one
  would fix it to the number of available CPUs.

- x:

  An object of class `"bootdist"` or `"density.bootdist"`.

- object:

  An object of class `"bootdist"`.

- main:

  an overall title for the plot: see
  [`title`](https://rdrr.io/r/graphics/title.html), default to
  `"Bootstrapped values of parameters"`.

- enhance:

  a logical to get an enhanced plot.

- trueval:

  when relevant, a numeric vector with the true value of parameters (for
  backfitting purposes).

- rampcol:

  colors to interpolate; must be a valid argument to
  [`colorRampPalette()`](https://rdrr.io/r/grDevices/colorRamp.html).

- nbgrid:

  Number of grid points in each direction. Can be scalar or a length-2
  integer vector.

- nbcol:

  An integer argument, the required number of colors

- ...:

  Further arguments to be passed to generic methods or `"bootdist"`
  objects for `density`.

- bw, adjust, kernel:

  resp. the smoothing bandwidth, the scaling factor, the kernel used,
  see [`density`](https://rdrr.io/r/stats/density.html).

- mar:

  A numerical vector of the form `c(bottom, left, top, right)`, see
  [`par`](https://rdrr.io/r/graphics/par.html).

- lty, col, lwd:

  resp. the line type, the color, the line width, see
  [`par`](https://rdrr.io/r/graphics/par.html).

## Details

Samples are drawn by parametric bootstrap (resampling from the
distribution fitted by
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md))
or nonparametric bootstrap (resampling with replacement from the data
set). On each bootstrap sample the function
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)
(or
[`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md),
[`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md),
[`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md)
according to the component `f$method` of the object of class
`"fitdist"`) is used to estimate bootstrapped values of parameters. When
that function fails to converge, `NA` values are returned. Medians and
2.5 and 97.5 percentiles are computed by removing `NA` values. The
medians and the 95 percent confidence intervals of parameters (2.5 and
97.5 percentiles) are printed in the summary. If inferior to the whole
number of iterations, the number of iterations for which the function
converges is also printed in the summary.

By default (when `enhance=FALSE`), the plot of an object of class
`"bootdist"` consists in a scatterplot or a matrix of scatterplots of
the bootstrapped values of parameters. It uses the function
[`stripchart`](https://rdrr.io/r/graphics/stripchart.html) when the
fitted distribution is characterized by only one parameter, the function
[`plot`](https://rdrr.io/r/graphics/plot.default.html) when there are
two paramters and the function
[`pairs`](https://rdrr.io/r/graphics/pairs.html) in other cases. In
these last cases, it provides a representation of the joint uncertainty
distribution of the fitted parameters.

When `enhance=TRUE`, a personalized plot version of
[`pairs`](https://rdrr.io/r/graphics/pairs.html) is used where upper
graphs are scatterplots and lower graphs are heatmap image using
[`image`](https://rdrr.io/r/graphics/image.html) based on a kernel based
estimator for the 2D density function (using `kde2d` from MASS package).
Arguments `rampcol`, `nbgrid`, `nbcol` can be used to customize the
plots. Defautls values are
`rampcol=c("green", "yellow", "orange", "red")`, `nbcol=100` (see
[`colorRampPalette()`](https://rdrr.io/r/grDevices/colorRamp.html)),
`nbgrid=100` (see `kde2d`). In addition, when fitting parameters on
simulated datasets for backtesting purposes, an additional argument
`trueval` can be used to plot a cross at the true value.

It is possible to accelerate the bootstrap using parallelization. We
recommend you to use `parallel = "multicore"`, or `parallel = "snow"` if
you work on Windows, and to fix `ncpus` to the number of available
processors.

`density` computes the empirical density of `bootdist` objects using the
[`density`](https://rdrr.io/r/stats/density.html) function (with
Gaussian kernel by default). It returns an object of class
`density.bootdist` for which `print` and `plot` methods are provided.

## Value

`bootdist` returns an object of class `"bootdist"`, a list with 6
components,

- estim:

  a data frame containing the bootstrapped values of parameters.

- converg:

  a vector containing the codes for convergence obtained if an iterative
  method is used to estimate parameters on each bootstraped data set
  (and 0 if a closed formula is used).

- method:

  A character string coding for the type of resampling : `"param"` for a
  parametric resampling and `"nonparam"` for a nonparametric resampling.

- nbboot:

  The number of samples drawn by bootstrap.

- CI:

  bootstrap medians and 95 percent confidence percentile intervals of
  parameters.

- fitpart:

  The object of class `"fitdist"` on which the bootstrap procedure was
  applied.

Generic functions:

- `print`:

  The print of a `"bootdist"` object shows the bootstrap parameter
  estimates. If inferior to the whole number of bootstrap iterations,
  the number of iterations for which the estimation converges is also
  printed.

- `summary`:

  The summary provides the median and 2.5 and 97.5 percentiles of each
  parameter. If inferior to the whole number of bootstrap iterations,
  the number of iterations for which the estimation converges is also
  printed in the summary.

- `plot`:

  The plot shows the bootstrap estimates with
  [`stripchart`](https://rdrr.io/r/graphics/stripchart.html) function
  for univariate parameters and
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) function for
  multivariate parameters.

- `density`:

  The density computes empirical densities and return an object of class
  `density.bootdist`.

## See also

See
[`fitdistrplus`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistrplus.md)
for an overview of the package.
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md),
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
[`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md),
[`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md),
[`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md),
[`quantile.bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/quantile.md)
for another generic function to calculate quantiles from the fitted
distribution and its bootstrap results and
[`CIcdfplot`](https://lbbe-software.github.io/fitdistrplus/reference/CIcdfplot.md)
for adding confidence intervals on quantiles to a CDF plot of the fitted
distribution.

Please visit the [Frequently Asked
Questions](https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html).

## References

Cullen AC and Frey HC (1999), *Probabilistic techniques in exposure
assessment*. Plenum Press, USA, pp. 181-241.

Delignette-Muller ML and Dutang C (2015), *fitdistrplus: An R Package
for Fitting Distributions*. Journal of Statistical Software, 64(4),
1-34, [doi:10.18637/jss.v064.i04](https://doi.org/10.18637/jss.v064.i04)
.

## Author

Marie-Laure Delignette-Muller and Christophe Dutang.

## Examples

``` r
# We choose a low number of bootstrap replicates in order to satisfy CRAN running times
# constraint.
# For practical applications, we recommend to use at least niter=501 or niter=1001.

set.seed(123) # here just to make random sampling reproducible

# (1) Fit of a gamma distribution to serving size data
# using default method (maximum likelihood estimation)
# followed by parametric bootstrap
#
data(groundbeef)
x1 <- groundbeef$serving
f1 <- fitdist(x1, "gamma")
b1 <- bootdist(f1, niter=51)
print(b1)
#> Parameter values obtained with parametric bootstrap 
#>      shape       rate
#> 1 4.212554 0.05996094
#> 2 4.340921 0.06000622
#> 3 4.073375 0.05552513
#> 4 3.845651 0.05427429
#> 5 4.194284 0.05601680
#> 6 3.812139 0.05111551
plot(b1)

plot(b1, enhance=TRUE)

summary(b1)
#> Parametric bootstrap medians and 95% percentile CI 
#>           Median       2.5%      97.5%
#> shape 4.15104594 3.34046628 4.90950540
#> rate  0.05679841 0.04571069 0.06890783
quantile(b1)
#> (original) estimated quantiles for each specified probability (non-censored data)
#>             p=0.1    p=0.2    p=0.3    p=0.4    p=0.5    p=0.6    p=0.7
#> estimate 32.16733 42.32692 50.91831 59.15298 67.62801 76.88308 87.67764
#>             p=0.8    p=0.9
#> estimate 101.5208 122.9543
#> Median of bootstrap estimates
#>             p=0.1    p=0.2    p=0.3    p=0.4    p=0.5    p=0.6    p=0.7
#> estimate 32.08309 42.39622 50.90926 59.00348 67.09501 76.57024 87.10971
#>             p=0.8    p=0.9
#> estimate 100.7966 121.6413
#> 
#> two-sided 95 % CI of each quantile
#>           p=0.1    p=0.2    p=0.3    p=0.4    p=0.5    p=0.6    p=0.7     p=0.8
#> 2.5 %  28.38291 38.60850 47.34299 55.47748 63.76341 72.42823 82.07640  94.18397
#> 97.5 % 35.50722 45.92833 54.68875 63.05277 71.83915 81.23027 92.03456 106.37599
#>           p=0.9
#> 2.5 %  113.3739
#> 97.5 % 129.3437
CIcdfplot(b1, CI.output = "quantile")

density(b1)
#> 
#> Bootstrap values for: gamma for 1 object(s) with 51 bootstrap values (original sample size 254).
plot(density(b1))


# (2) non parametric bootstrap on the same fit
#
b1b <- bootdist(f1, bootmethod="nonparam", niter=51) 
summary(b1b)
#> Nonparametric bootstrap medians and 95% percentile CI 
#>         Median       2.5%      97.5%
#> shape 3.956460 3.50471143 4.81001129
#> rate  0.054112 0.04729992 0.06716172
quantile(b1b)
#> (original) estimated quantiles for each specified probability (non-censored data)
#>             p=0.1    p=0.2    p=0.3    p=0.4    p=0.5    p=0.6    p=0.7
#> estimate 32.16733 42.32692 50.91831 59.15298 67.62801 76.88308 87.67764
#>             p=0.8    p=0.9
#> estimate 101.5208 122.9543
#> Median of bootstrap estimates
#>             p=0.1    p=0.2    p=0.3    p=0.4    p=0.5    p=0.6    p=0.7
#> estimate 31.80216 41.79758 50.50144 58.77236 67.23737 76.63013 87.58129
#>             p=0.8    p=0.9
#> estimate 101.2415 122.7574
#> 
#> two-sided 95 % CI of each quantile
#>           p=0.1    p=0.2    p=0.3    p=0.4    p=0.5    p=0.6    p=0.7     p=0.8
#> 2.5 %  29.20585 39.00841 47.27673 55.44208 63.53759 72.05392 81.78703  94.34323
#> 97.5 % 35.82375 46.10281 54.74061 62.95224 71.34868 80.56793 91.36288 105.53125
#>           p=0.9
#> 2.5 %  113.7319
#> 97.5 % 128.7492


# (3) Fit of a normal distribution on acute toxicity values of endosulfan in log10 for
# nonarthropod invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5 percent quantile value of 
# the fitted distribution, what is called the 5 percent hazardous concentration (HC5)
# in ecotoxicology, with its two-sided 95 percent confidence interval calculated by 
# parametric bootstrap
#
data(endosulfan)
ATV <- subset(endosulfan, group == "NonArthroInvert")$ATV
log10ATV <- log10(subset(endosulfan, group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV, "norm")
bln <- bootdist(fln, bootmethod = "param", niter=51)
quantile(bln, probs = c(0.05, 0.1, 0.2))
#> (original) estimated quantiles for each specified probability (non-censored data)
#>            p=0.05    p=0.1  p=0.2
#> estimate 1.744227 2.080093 2.4868
#> Median of bootstrap estimates
#>           p=0.05    p=0.1    p=0.2
#> estimate 1.86876 2.180432 2.615232
#> 
#> two-sided 95 % CI of each quantile
#>         p=0.05    p=0.1    p=0.2
#> 2.5 %  1.37566 1.729433 2.186051
#> 97.5 % 2.79308 2.986272 3.237189

# (4) comparison of sequential and parallel versions of bootstrap
# to be tried with a greater number of iterations (1001 or more)
#
# \donttest{
niter <- 1001
data(groundbeef)
x1 <- groundbeef$serving
f1 <- fitdist(x1, "gamma")

# sequential version
ptm <- proc.time()
summary(bootdist(f1, niter = niter))
#> Parametric bootstrap medians and 95% percentile CI 
#>           Median       2.5%      97.5%
#> shape 4.02759987 3.38936853 4.79038017
#> rate  0.05465955 0.04584153 0.06559032
proc.time() - ptm
#>    user  system elapsed 
#>   4.166   0.000   4.166 

# parallel version using snow
require("parallel")
#> Loading required package: parallel
ptm <- proc.time()
summary(bootdist(f1, niter = niter, parallel = "snow", ncpus = 2))
#> Parametric bootstrap medians and 95% percentile CI 
#>           Median       2.5%      97.5%
#> shape 4.03717027 3.42095789 4.78943116
#> rate  0.05477864 0.04620101 0.06603422
proc.time() - ptm
#>    user  system elapsed 
#>   0.038   0.002   3.885 

# parallel version using multicore (not available on Windows)
ptm <- proc.time()
summary(bootdist(f1, niter = niter, parallel = "multicore", ncpus = 2))
#> Parametric bootstrap medians and 95% percentile CI 
#>           Median       2.5%     97.5%
#> shape 4.02586191 3.45250533 4.7634281
#> rate  0.05466294 0.04663766 0.0655417
proc.time() - ptm
#>    user  system elapsed 
#>   4.262   0.258   2.285 
# }
```
