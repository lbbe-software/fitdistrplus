# Matching moment fit of univariate distributions

Fit of univariate distributions by matching moments (raw or centered)
for non censored data.

## Usage

``` r
mmedist(data, distr, order, memp, start = NULL, fix.arg = NULL, optim.method = "default", 
  lower = -Inf, upper = Inf, custom.optim = NULL, weights = NULL, silent = TRUE, 
  gradient = NULL, checkstartfix=FALSE, calcvcov=FALSE, ...)
```

## Arguments

- data:

  A numeric vector for non censored data.

- distr:

  A character string `"name"` naming a distribution (see 'details').

- order:

  A numeric vector for the moment order(s). The length of this vector
  must be equal to the number of parameters to estimate.

- memp:

  A function implementing empirical moments, raw or centered but has to
  be consistent with `distr` argument (and `weights` argument). See
  details below.

- start:

  A named list giving the initial values of parameters of the named
  distribution or a function of data computing initial values and
  returning a named list. This argument may be omitted (default) for
  some distributions for which reasonable starting values are computed
  (see the 'details' section of
  [`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)).

- fix.arg:

  An optional named list giving the values of fixed parameters of the
  named distribution or a function of data computing (fixed) parameter
  values and returning a named list. Parameters with fixed value are
  thus NOT estimated.

- optim.method:

  `"default"` or optimization method to pass to
  [`optim`](https://rdrr.io/r/stats/optim.html).

- lower:

  Left bounds on the parameters for the `"L-BFGS-B"` method (see
  [`optim`](https://rdrr.io/r/stats/optim.html)).

- upper:

  Right bounds on the parameters for the `"L-BFGS-B"` method (see
  [`optim`](https://rdrr.io/r/stats/optim.html)).

- custom.optim:

  a function carrying the optimization .

- weights:

  an optional vector of weights to be used in the fitting process.
  Should be `NULL` or a numeric vector with strictly positive integers
  (typically the number of occurences of each observation). If
  non-`NULL`, weighted MME is used, otherwise ordinary MME.

- silent:

  A logical to remove or show warnings when bootstraping.

- gradient:

  A function to return the gradient of the squared difference for the
  `"BFGS"`, `"CG"` and `"L-BFGS-B"` methods. If it is `NULL`, a
  finite-difference approximation will be used, see details.

- checkstartfix:

  A logical to test starting and fixed values. Do not change it.

- calcvcov:

  A logical indicating if (asymptotic) covariance matrix is required.

- ...:

  further arguments passed to the
  [`optim`](https://rdrr.io/r/stats/optim.html),
  [`constrOptim`](https://rdrr.io/r/stats/constrOptim.html) or
  `custom.optim` function.

## Details

The argument `distr` can be one of the base R distributions: `"norm"`,
`"lnorm"`, `"exp"` and `"pois"`, `"gamma"`, `"logis"`, `"nbinom"` ,
`"geom"`, `"beta"` and `"unif"`. In that case, no other arguments than
`data` and `distr` are required, because the estimate is computed by a
closed-form formula. For distributions characterized by one parameter
(`"geom"`, `"pois"` and `"exp"`), this parameter is simply estimated by
matching theoretical and observed means, and for distributions
characterized by two parameters, these parameters are estimated by
matching theoretical and observed means and variances (Vose, 2000). Note
that for these closed-form formula, `fix.arg` cannot be used and `start`
is ignored.

The argument `distr` can also be the distribution name as long as a
corresponding `mdistr` function exists, e.g. `"pareto"` if `"mpareto"`
exists. In that case arguments arguments `order` and `memp` have to be
supplied in order to carry out the matching numerically, by minimization
of the sum of squared differences between observed and theoretical
moments. Optionnally other arguments can be supplied to control
optimization (see the 'details' section of
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)
for details about arguments for the control of optimization). In that
case, `fix.arg` can be used and `start` is taken into account.

For non closed-form estimators, `memp` must be provided to compute
empirical moments. When `weights=NULL`, this function must have two
arguments `x, order`: `x` the numeric vector of the data and `order` the
order of the moment. When `weights!=NULL`, this function must have three
arguments `x, order, weights`: `x` the numeric vector of the data,
`order` the order of the moment, `weights` the numeric vector of
weights. See examples below.

Optionally, a vector of `weights` can be used in the fitting process. By
default (when `weigths=NULL`), ordinary MME is carried out, otherwise
the specified weights are used to compute (raw or centered) weighted
moments. For closed-form estimators, weighted mean and variance are
computed by `wtdmean` and `wtdvar` from the `Hmisc` package. When a
numerical minimization is used, weighted are expected to be computed by
the `memp` function. It is not yet possible to take into account
weighths in functions `plotdist`, `plotdistcens`, `plot.fitdist`,
`plot.fitdistcens`, `cdfcomp`, `cdfcompcens`, `denscomp`, `ppcomp`,
`qqcomp`, `gofstat` and `descdist` (developments planned in the future).

This function is not intended to be called directly but is internally
called in
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
and
[`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md)
when used with the matching moments method.

Since Version 1.2-0, `mmedist` automatically computes the asymptotic
covariance matrix using I. Ibragimov and R. Has'minskii (1981), hence
the theoretical moments `mdist` should be defined up to an order which
equals to twice the maximal order given `order`. For instance, the
normal distribution, we fit against the expectation and the variance and
we need to have `mnorm` up to order \\2\times2=4\\.

## Value

`mmedist` returns a list with following components,

- estimate:

  the parameter estimates.

- convergence:

  an integer code for the convergence of
  [`optim`](https://rdrr.io/r/stats/optim.html) defined as below or
  defined by the user in the user-supplied optimization function. `0`
  indicates successful convergence. `1` indicates that the iteration
  limit of [`optim`](https://rdrr.io/r/stats/optim.html) has been
  reached. `10` indicates degeneracy of the Nealder-Mead simplex. `100`
  indicates that [`optim`](https://rdrr.io/r/stats/optim.html)
  encountered an internal error.

- value:

  the minimal value reached for the criterion to minimize.

- hessian:

  a symmetric matrix computed by
  [`optim`](https://rdrr.io/r/stats/optim.html) as an estimate of the
  Hessian at the solution found or computed in the user-supplied
  optimization function.

- optim.function:

  (if appropriate) the name of the optimization function used for
  maximum likelihood.

- optim.method:

  (if appropriate) when [`optim`](https://rdrr.io/r/stats/optim.html) is
  used, the name of the algorithm used, the field `method` of the
  `custom.optim` function otherwise.

- fix.arg:

  the named list giving the values of parameters of the named
  distribution that must kept fixed rather than estimated by maximum
  likelihood or `NULL` if there are no such parameters.

- fix.arg.fun:

  the function used to set the value of `fix.arg` or `NULL`.

- weights:

  the vector of weigths used in the estimation process or `NULL`.

- counts:

  A two-element integer vector giving the number of calls to the
  log-likelihood function and its gradient respectively. This excludes
  those calls needed to compute the Hessian, if requested, and any calls
  to log-likelihood function to compute a finite-difference
  approximation to the gradient. `counts` is returned by
  [`optim`](https://rdrr.io/r/stats/optim.html) or the user-supplied
  function or set to `NULL`.

- optim.message:

  A character string giving any additional information returned by the
  optimizer, or `NULL`. To understand exactly the message, see the
  source code.

- loglik:

  the log-likelihood value.

- method:

  either `"closed formula"` or the name of the optimization method.

- order:

  the order of the moment(s) matched.

- memp:

  the empirical moment function.

## See also

See
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
[`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md),
[`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md),
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md),[`fitdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.md),
[`optim`](https://rdrr.io/r/stats/optim.html),
[`bootdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/bootdistcens.md)
and
[`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md).

Please visit the [Frequently Asked
Questions](https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html).

## References

I. Ibragimov and R. Has'minskii (1981), *Statistical Estimation -
Asymptotic Theory*, Springer-Verlag,
[doi:10.1007/978-1-4899-0027-2](https://doi.org/10.1007/978-1-4899-0027-2)

Evans M, Hastings N and Peacock B (2000), *Statistical distributions*.
John Wiley and Sons Inc,
[doi:10.1002/9780470627242](https://doi.org/10.1002/9780470627242) .

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
set.seed(123) # here just to make random sampling reproducible

# (1) basic fit of a normal distribution with moment matching estimation
#

n <- 100
x1 <- rnorm(n=n)
mmedist(x1, "norm")
#> $estimate
#>       mean         sd 
#> 0.09040591 0.90824033 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> NULL
#> 
#> $hessian
#> NULL
#> 
#> $optim.function
#> NULL
#> 
#> $opt.meth
#> NULL
#> 
#> $fix.arg
#> NULL
#> 
#> $fix.arg.fun
#> NULL
#> 
#> $weights
#> NULL
#> 
#> $counts
#> NULL
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -132.2692
#> 
#> $method
#> [1] "closed formula"
#> 
#> $order
#> [1] 1 2
#> 
#> $memp
#> NULL
#> 
#> $vcov
#> NULL
#> 

#weighted
w <- c(rep(1, n/2), rep(10, n/2))
mmedist(x1, "norm", weights=w)$estimate
#> Warning: weights are not taken into account in the default initial values
#>      mean        sd 
#> 0.1362260 0.8995989 


# (2) fit a discrete distribution (Poisson)
#

x2 <- rpois(n=30,lambda = 2)
mmedist(x2, "pois")
#> $estimate
#>   lambda 
#> 2.066667 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> NULL
#> 
#> $hessian
#> NULL
#> 
#> $optim.function
#> NULL
#> 
#> $opt.meth
#> NULL
#> 
#> $fix.arg
#> NULL
#> 
#> $fix.arg.fun
#> NULL
#> 
#> $weights
#> NULL
#> 
#> $counts
#> NULL
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -49.91188
#> 
#> $method
#> [1] "closed formula"
#> 
#> $order
#> [1] 1
#> 
#> $memp
#> NULL
#> 
#> $vcov
#> NULL
#> 

# (3) fit a finite-support distribution (beta)
#

x3 <- rbeta(n=100,shape1=5, shape2=10)
mmedist(x3, "beta")
#> $estimate
#>    shape1    shape2 
#>  5.157206 10.904341 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> NULL
#> 
#> $hessian
#> NULL
#> 
#> $optim.function
#> NULL
#> 
#> $opt.meth
#> NULL
#> 
#> $fix.arg
#> NULL
#> 
#> $fix.arg.fun
#> NULL
#> 
#> $weights
#> NULL
#> 
#> $counts
#> NULL
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] 79.14058
#> 
#> $method
#> [1] "closed formula"
#> 
#> $order
#> [1] 1 2
#> 
#> $memp
#> NULL
#> 
#> $vcov
#> NULL
#> 


# (4) fit a Pareto distribution
#

# \donttest{
  require("actuar")
  #simulate a sample
  x4  <-  rpareto(1000, 6, 2)

  #empirical raw moment
  memp  <-  function(x, order) mean(x^order)
  memp2 <- function(x, order, weights) sum(x^order * weights)/sum(weights)

  #fit by MME
  mmedist(x4, "pareto", order=c(1, 2), memp=memp, 
    start=list(shape=10, scale=10), lower=1, upper=Inf)
#> $estimate
#>    shape    scale 
#> 5.818372 1.998979 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 1.590384e-13
#> 
#> $hessian
#> NULL
#> 
#> $optim.function
#> [1] "constrOptim"
#> 
#> $optim.method
#> [1] "Nelder-Mead"
#> 
#> $fix.arg
#> NULL
#> 
#> $fix.arg.fun
#> NULL
#> 
#> $weights
#> NULL
#> 
#> $counts
#> function gradient 
#>      911       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -102.229
#> 
#> $method
#> [1] "default"
#> 
#> $order
#> [1] 1 2
#> 
#> $memp
#> function (x, order) 
#> mean(x^order)
#> <environment: 0x55e416b5fdd8>
#> 
#> $vcov
#> NULL
#> 
  #fit by weighted MME
  w <- rep(1, length(x4))
  w[x4 < 1] <- 2
  mmedist(x4, "pareto", order=c(1, 2), memp=memp2, weights=w,
    start=list(shape=10, scale=10), lower=1, upper=Inf)
#> Warning: weights are not taken into account in the default initial values
#> $estimate
#>    shape    scale 
#> 7.482125 2.284692 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 9.52255e-13
#> 
#> $hessian
#> NULL
#> 
#> $optim.function
#> [1] "constrOptim"
#> 
#> $optim.method
#> [1] "Nelder-Mead"
#> 
#> $fix.arg
#> NULL
#> 
#> $fix.arg.fun
#> NULL
#> 
#> $weights
#>    [1] 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2
#>   [38] 2 2 1 2 2 1 2 2 1 2 2 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 1 2 2 2 2 1 2 2
#>   [75] 2 1 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 1 1 2 2
#>  [112] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1
#>  [149] 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 1 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2
#>  [186] 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 1 2 1 2
#>  [223] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2
#>  [260] 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 2 2 1 2 2 2 2 1
#>  [297] 2 2 2 2 2 1 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2
#>  [334] 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2
#>  [371] 2 2 2 1 2 1 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 1
#>  [408] 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#>  [445] 1 2 1 2 2 2 2 2 2 2 2 2 2 1 1 1 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2
#>  [482] 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2
#>  [519] 2 1 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#>  [556] 2 1 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#>  [593] 2 1 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2
#>  [630] 2 2 2 2 2 2 1 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2
#>  [667] 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2
#>  [704] 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2
#>  [741] 2 2 1 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 1 2 2 2
#>  [778] 2 1 2 2 2 2 1 2 2 2 2 1 2 2 2 2 2 2 1 2 2 2 2 2 1 2 2 1 2 2 2 2 2 1 2 2 1
#>  [815] 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2
#>  [852] 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2
#>  [889] 2 2 2 2 2 2 2 2 2 2 1 1 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#>  [926] 2 1 2 2 2 2 1 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 2 2 2 1 2 2 2 2 2 1 2 2 2 2 2
#>  [963] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 1 2 2 2 2
#> [1000] 2
#> 
#> $counts
#> function gradient 
#>     1560       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] 95.23506
#> 
#> $method
#> [1] "default"
#> 
#> $order
#> [1] 1 2
#> 
#> $memp
#> function (x, order, weights) 
#> sum(x^order * weights)/sum(weights)
#> <environment: 0x55e416b5fdd8>
#> 
#> $vcov
#> NULL
#> 
# }
```
