# Bootstrap simulation of uncertainty for censored data

Uses nonparametric bootstrap resampling in order to simulate uncertainty
in the parameters of the distribution fitted to censored data.

## Usage

``` r
bootdistcens(f, niter = 1001, silent = TRUE, 
      parallel = c("no", "snow", "multicore"), ncpus)
# S3 method for class 'bootdistcens'
print(x, ...)
# S3 method for class 'bootdistcens'
plot(x, ...)
# S3 method for class 'bootdistcens'
summary(object, ...)
# S3 method for class 'bootdistcens'
density(..., bw = nrd0, adjust = 1, kernel = "gaussian")
# S3 method for class 'density.bootdistcens'
plot(x, mar=c(4,4,2,1), lty=NULL, col=NULL, lwd=NULL, ...)
# S3 method for class 'density.bootdistcens'
print(x, ...)
```

## Arguments

- f:

  An object of class `"fitdistcens"`, output of the
  [`fitdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.md)
  function.

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

  An object of class `"bootdistcens"`.

- object:

  An object of class `"bootdistcens"`.

- ...:

  Further arguments to be passed to generic methods or `"bootdistcens"`
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

Samples are drawn by nonparametric bootstrap (resampling with
replacement from the data set). On each bootstrap sample the function
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)
is used to estimate bootstrapped values of parameters. When
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)
fails to converge, `NA` values are returned. Medians and 2.5 and 97.5
percentiles are computed by removing `NA` values. The medians and the 95
percent confidence intervals of parameters (2.5 and 97.5 percentiles)
are printed in the summary. If inferior to the whole number of
iterations, the number of iterations for which
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)
converges is also printed in the summary.

The plot of an object of class `"bootdistcens"` consists in a
scatterplot or a matrix of scatterplots of the bootstrapped values of
parameters. It uses the function
[`stripchart`](https://rdrr.io/r/graphics/stripchart.html) when the
fitted distribution is characterized by only one parameter, and the
function [`plot`](https://rdrr.io/r/graphics/plot.default.html) in other
cases. In these last cases, it provides a representation of the joint
uncertainty distribution of the fitted parameters.

It is possible to accelerate the bootstrap using parallelization. We
recommend you to use `parallel = "multicore"`, or `parallel = "snow"` if
you work on Windows, and to fix `ncpus` to the number of available
processors.

`density` computes the empirical density of `bootdistcens` objects using
the [`density`](https://rdrr.io/r/stats/density.html) function (with
Gaussian kernel by default). It returns an object of class
`density.bootdistcens` for which `print` and `plot` methods are
provided.

## Value

`bootdistcens` returns an object of class `"bootdistcens"`, a list with
6 components,

- estim:

  a data frame containing the bootstrapped values of parameters.

- converg:

  a vector containing the codes for convergence of the iterative method
  used to estimate parameters on each bootstraped data set.

- method:

  A character string coding for the type of resampling : in this case
  `"nonparam"` as it is the only available method for censored data.

- nbboot:

  The number of samples drawn by bootstrap.

- CI:

  bootstrap medians and 95 percent confidence percentile intervals of
  parameters.

- fitpart:

  The object of class `"fitdistcens"` on which the bootstrap procedure
  was applied.

Generic functions:

- `print`:

  The print of a `"bootdistcens"` object shows the bootstrap parameter
  estimates. If inferior to the whole number of bootstrap iterations,
  the number of iterations for which the estimation converges is also
  printed.

- `summary`:

  The summary provides the median and 2.5 and 97.5 percentiles of each
  parameter. If inferior to the whole number of bootstrap iterations,
  the number of iterations for which the estimation converges is also
  printed in the summary.

- `plot`:

  The plot shows the bootstrap estimates with the
  [`stripchart`](https://rdrr.io/r/graphics/stripchart.html) function
  for univariate parameters and
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) function for
  multivariate parameters.

- `density`:

  The density computes empirical densities and return an object of class
  `density.bootdistcens`.

## See also

See
[`fitdistrplus`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistrplus.md)
for an overview of the package.
[`fitdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.md),
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
[`quantile.bootdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/quantile.md)
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

# (1) Fit of a normal distribution to fluazinam data in log10
# followed by nonparametric bootstrap and calculation of quantiles
# with 95 percent confidence intervals
#
data(fluazinam)
(d1 <-log10(fluazinam))
#>         left     right
#> 1  0.5797836 0.5797836
#> 2  1.5263393 1.5263393
#> 3  1.9395193 1.9395193
#> 4  3.2304489        NA
#> 5  2.8061800 2.8061800
#> 6  3.0625820        NA
#> 7  2.0530784 2.0530784
#> 8  2.1105897 2.1105897
#> 9  2.7678976 2.7678976
#> 10 3.2685780        NA
#> 11 0.2041200 0.2041200
#> 12 0.6812412 0.6812412
#> 13 1.9138139 1.9138139
#> 14 2.1903317 2.1903317
f1 <- fitdistcens(d1, "norm")
b1 <- bootdistcens(f1, niter = 51)
b1
#> Parameter values obtained with nonparametric bootstrap 
#>        mean        sd
#> 1  2.602889 1.4296898
#> 2  2.306756 0.9805219
#> 3  2.123152 1.0523904
#> 4  2.104933 1.1213470
#> 5  2.189861 1.2861038
#> 6  2.394168 1.4639464
#> 7  2.490409 0.7827611
#> 8  2.349840 0.5102881
#> 9  2.672415 1.6120354
#> 10 2.034794 0.8691045
#> 11 2.122176 1.2861032
#> 12 2.400458 0.6773500
#> 13 2.406236 1.1745370
#> 14 2.506662 1.6512832
#> 15 2.182508 1.0427324
#> 16 2.220548 0.6878484
#> 17 2.295152 1.3031345
#> 18 2.688828 1.4270509
#> 19 2.413936 1.0708897
#> 20 2.581404 1.6380040
#> 21 2.289350 1.0410812
#> 22 2.522076 0.8808777
#> 23 2.007564 0.9032233
#> 24 1.569873 0.7254079
#> 25 2.662200 1.7368653
#> 26 2.392542 0.7067813
#> 27 1.419076 0.8639953
#> 28 1.921136 0.8581193
#> 29 3.286352 1.1160621
#> 30 1.754207 1.0323776
#> 31 2.100860 1.1962162
#> 32 1.817393 1.0456629
#> 33 2.053392 1.4666924
#> 34 2.222665 1.2763759
#> 35 1.695874 1.1357282
#> 36 2.869931 1.2798787
#> 37 2.763266 1.1941164
#> 38 2.041762 1.0814692
#> 39 2.461076 0.8221329
#> 40 2.084426 1.0984889
#> 41 2.105346 1.4698270
#> 42 2.128899 1.1328278
#> 43 2.368866 1.5131334
#> 44 1.985509 0.9494199
#> 45 1.935726 0.9611150
#> 46 2.020270 1.1862202
#> 47 2.430016 1.4696352
#> 48 2.413309 1.1079596
#> 49 2.161517 1.0263700
#> 50 1.980760 0.9794176
#> 51 2.061176 1.0257701
summary(b1)
#> Nonparametric bootstrap medians and 95% percentile CI 
#>        Median      2.5%    97.5%
#> mean 2.220548 1.6013733 2.843264
#> sd   1.098489 0.6799746 1.647963
plot(b1)

quantile(b1)
#> (original) estimated quantiles for each specified probability (censored data)
#>              p=0.1    p=0.2    p=0.3   p=0.4    p=0.5    p=0.6    p=0.7
#> estimate 0.6655064 1.179033 1.549321 1.86572 2.161449 2.457179 2.773577
#>             p=0.8    p=0.9
#> estimate 3.143865 3.657392
#> Median of bootstrap estimates
#>              p=0.1    p=0.2    p=0.3    p=0.4    p=0.5   p=0.6    p=0.7
#> estimate 0.7255863 1.198924 1.611787 1.918335 2.220548 2.47912 2.796609
#>             p=0.8    p=0.9
#> estimate 3.131984 3.563346
#> 
#> two-sided 95 % CI of each quantile
#>           p=0.1     p=0.2    p=0.3    p=0.4    p=0.5    p=0.6    p=0.7    p=0.8
#> 2.5 %  0.226360 0.7597641 1.122590 1.391605 1.601373 1.811142 2.035571 2.291063
#> 97.5 % 1.655009 1.8981832 2.183339 2.524443 2.843264 3.171196 3.565035 4.100272
#>           p=0.9
#> 2.5 %  2.645698
#> 97.5 % 4.732902
CIcdfplot(b1, CI.output = "quantile")

plot(density(b1))
#> List of 1
#>  $ :List of 6
#>   ..$ estim  :'data.frame':  51 obs. of  2 variables:
#>   .. ..$ mean: num [1:51] 2.6 2.31 2.12 2.1 2.19 ...
#>   .. ..$ sd  : num [1:51] 1.43 0.981 1.052 1.121 1.286 ...
#>   ..$ converg: num [1:51] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ method : chr "nonparam"
#>   ..$ nbboot : num 51
#>   ..$ CI     : num [1:2, 1:3] 2.22 1.1 1.6 0.68 2.84 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:2] "mean" "sd"
#>   .. .. ..$ : chr [1:3] "Median" "2.5%" "97.5%"
#>   ..$ fitpart:List of 17
#>   .. ..$ estimate   : Named num [1:2] 2.16 1.17
#>   .. .. ..- attr(*, "names")= chr [1:2] "mean" "sd"
#>   .. ..$ method     : chr "mle"
#>   .. ..$ sd         : Named num [1:2] 0.322 0.263
#>   .. .. ..- attr(*, "names")= chr [1:2] "mean" "sd"
#>   .. ..$ cor        : num [1:2, 1:2] 1 0.135 0.135 1
#>   .. .. ..- attr(*, "dimnames")=List of 2
#>   .. .. .. ..$ : chr [1:2] "mean" "sd"
#>   .. .. .. ..$ : chr [1:2] "mean" "sd"
#>   .. ..$ vcov       : num [1:2, 1:2] 0.1039 0.0114 0.0114 0.0692
#>   .. .. ..- attr(*, "dimnames")=List of 2
#>   .. .. .. ..$ : chr [1:2] "mean" "sd"
#>   .. .. .. ..$ : chr [1:2] "mean" "sd"
#>   .. ..$ loglik     : num -20.4
#>   .. ..$ aic        : num 44.8
#>   .. ..$ bic        : num 46.1
#>   .. ..$ n          : int 14
#>   .. ..$ censdata   :'data.frame':   14 obs. of  2 variables:
#>   .. .. ..$ left : num [1:14] 0.58 1.53 1.94 3.23 2.81 ...
#>   .. .. ..$ right: num [1:14] 0.58 1.53 1.94 NA 2.81 ...
#>   .. ..$ distname   : chr "norm"
#>   .. ..$ fix.arg    : NULL
#>   .. ..$ fix.arg.fun: NULL
#>   .. ..$ dots       : NULL
#>   .. ..$ convergence: int 0
#>   .. ..$ discrete   : logi FALSE
#>   .. ..$ weights    : NULL
#>   .. ..- attr(*, "class")= chr "fitdistcens"
#>   ..- attr(*, "class")= chr "bootdistcens"
#> NULL


# (2) Estimation of the mean of the normal distribution 
# by maximum likelihood with the standard deviation fixed at 1 
# using the argument fix.arg
# followed by nonparametric bootstrap 
# and calculation of quantiles with 95 percent confidence intervals
#
f1b <- fitdistcens(d1, "norm", start = list(mean = 1),fix.arg = list(sd = 1))
b1b <- bootdistcens(f1b, niter = 51)
summary(b1b)
#> Nonparametric bootstrap medians and 95% percentile CI 
#>   Median     2.5%    97.5% 
#> 2.131884 1.675051 2.718799 
plot(b1b)
quantile(b1b)
#> (original) estimated quantiles for each specified probability (censored data)
#>              p=0.1    p=0.2    p=0.3    p=0.4    p=0.5    p=0.6    p=0.7
#> estimate 0.8527461 1.292676 1.609897 1.880951 2.134298 2.387645 2.658698
#>             p=0.8    p=0.9
#> estimate 2.975919 3.415849
#> Median of bootstrap estimates
#>             p=0.1    p=0.2    p=0.3    p=0.4    p=0.5    p=0.6    p=0.7
#> estimate 0.850332 1.290262 1.607483 1.878536 2.131884 2.385231 2.656284
#>             p=0.8    p=0.9
#> estimate 2.973505 3.413435
#> 
#> two-sided 95 % CI of each quantile
#>            p=0.1    p=0.2    p=0.3    p=0.4    p=0.5    p=0.6    p=0.7    p=0.8
#> 2.5 %  0.3934997 0.833430 1.150651 1.421704 1.675051 1.928398 2.199452 2.516672
#> 97.5 % 1.4372475 1.877178 2.194399 2.465452 2.718799 2.972146 3.243200 3.560420
#>           p=0.9
#> 2.5 %  2.956603
#> 97.5 % 4.000351

# (3) comparison of sequential and parallel versions of bootstrap
# to be tried with a greater number of iterations (1001 or more)
#
# \donttest{
niter <- 1001
data(fluazinam)
d1 <-log10(fluazinam)
f1 <- fitdistcens(d1, "norm")

# sequential version
ptm <- proc.time()
summary(bootdistcens(f1, niter = niter))
#> Nonparametric bootstrap medians and 95% percentile CI 
#>        Median     2.5%    97.5%
#> mean 2.164973 1.554776 2.864004
#> sd   1.098434 0.686322 1.684879
proc.time() - ptm
#>    user  system elapsed 
#>   4.922   0.000   4.925 

# parallel version using snow
require("parallel")
ptm <- proc.time()
summary(bootdistcens(f1, niter = niter, parallel = "snow", ncpus = 2))
#> Nonparametric bootstrap medians and 95% percentile CI 
#>        Median      2.5%    97.5%
#> mean 2.173299 1.5614412 2.883032
#> sd   1.110030 0.6555083 1.690022
proc.time() - ptm
#>    user  system elapsed 
#>   0.008   0.001   3.551 

# parallel version using multicore (not available on Windows)
ptm <- proc.time()
summary(bootdistcens(f1, niter = niter, parallel = "multicore", ncpus = 2))
#> Nonparametric bootstrap medians and 95% percentile CI 
#>        Median      2.5%    97.5%
#> mean 2.174916 1.5660318 2.827636
#> sd   1.117121 0.7103254 1.667564
proc.time() - ptm
#>    user  system elapsed 
#>   5.098   0.277   2.704 
# }

```
