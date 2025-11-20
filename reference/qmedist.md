# Quantile matching fit of univariate distributions

Fit of univariate distribution by matching quantiles for non censored
data.

## Usage

``` r
qmedist(data, distr, probs, start = NULL, fix.arg = NULL, qtype = 7, 
    optim.method = "default", lower = -Inf, upper = Inf, 
    custom.optim = NULL, weights = NULL, silent = TRUE, gradient = NULL, 
    checkstartfix=FALSE, calcvcov=FALSE, ...)
```

## Arguments

- data:

  A numeric vector for non censored data.

- distr:

  A character string `"name"` naming a distribution for which the
  corresponding quantile function `qname` and the corresponding density
  distribution `dname` must be classically defined.

- probs:

  A numeric vector of the probabilities for which the quantile matching
  is done. The length of this vector must be equal to the number of
  parameters to estimate.

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

- qtype:

  The quantile type used by the R
  [`quantile`](https://lbbe-software.github.io/fitdistrplus/reference/quantile.md)
  function to compute the empirical quantiles, (default 7 corresponds to
  the default quantile method in R).

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

  a function carrying the optimization.

- weights:

  an optional vector of weights to be used in the fitting process.
  Should be `NULL` or a numeric vector with strictly positive integers
  (typically the number of occurences of each observation). If
  non-`NULL`, weighted QME is used, otherwise ordinary QME.

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
  (currently ignored)

- ...:

  further arguments passed to the
  [`optim`](https://rdrr.io/r/stats/optim.html),
  [`constrOptim`](https://rdrr.io/r/stats/constrOptim.html) or
  `custom.optim` function.

## Details

The `qmedist` function carries out the quantile matching numerically, by
minimization of the sum of squared differences between observed and
theoretical quantiles. Note that for discrete distribution, the sum of
squared differences is a step function and consequently, the optimum is
not unique, see the FAQ.

The optimization process is the same as
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
see the 'details' section of that function.

Optionally, a vector of `weights` can be used in the fitting process. By
default (when `weigths=NULL`), ordinary QME is carried out, otherwise
the specified weights are used to compute weighted quantiles used in the
squared differences. Weigthed quantiles are computed by `wtdquantile`
from the `Hmisc` package. It is not yet possible to take into account
weighths in functions `plotdist`, `plotdistcens`, `plot.fitdist`,
`plot.fitdistcens`, `cdfcomp`, `cdfcompcens`, `denscomp`, `ppcomp`,
`qqcomp`, `gofstat` and `descdist` (developments planned in the future).

This function is not intended to be called directly but is internally
called in
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
and
[`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md).

## Value

`qmedist` returns a list with following components,

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

  the name of the optimization function used for maximum likelihood.

- optim.method:

  when [`optim`](https://rdrr.io/r/stats/optim.html) is used, the name
  of the algorithm used, the field `method` of the `custom.optim`
  function otherwise.

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

- probs:

  the probability vector on which quantiles are matched.

## See also

See
[`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md),
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
[`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md),
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
for other estimation methods and
[`quantile`](https://lbbe-software.github.io/fitdistrplus/reference/quantile.md)
for empirical quantile estimation in R.

Please visit the [Frequently Asked
Questions](https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html).

## References

Klugman SA, Panjer HH and Willmot GE (2012), *Loss Models: From Data to
Decissions*, 4th edition. Wiley Series in Statistics for Finance,
Business and Economics, p. 253,
[doi:10.1198/tech.2006.s409](https://doi.org/10.1198/tech.2006.s409) .

Delignette-Muller ML and Dutang C (2015), *fitdistrplus: An R Package
for Fitting Distributions*. Journal of Statistical Software, 64(4),
1-34, [doi:10.18637/jss.v064.i04](https://doi.org/10.18637/jss.v064.i04)
.

## Author

Christophe Dutang and Marie Laure Delignette-Muller.

## Examples

``` r
# (1) basic fit of a normal distribution 
#

set.seed(1234)
x1 <- rnorm(n=100)
qmedist(x1, "norm", probs=c(1/3, 2/3))
#> $estimate
#>       mean         sd 
#> -0.3025734  0.8521385 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 2.427759e-10
#> 
#> $hessian
#>               mean            sd
#> mean  2.000000e+00 -2.784663e-14
#> sd   -2.784663e-14  3.710520e-01
#> 
#> $optim.function
#> [1] "optim"
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
#>       57       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -146.1278
#> 
#> $probs
#> [1] 0.3333333 0.6666667
#> 


# (2) defining your own distribution functions, here for the Gumbel 
# distribution for other distributions, see the CRAN task view dedicated 
# to probability distributions

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
qgumbel <- function(p, a, b) a - b*log(-log(p))
qmedist(x1, "gumbel", probs=c(1/3, 2/3), start=list(a=10,b=5))
#> Error in checkparamlist(arg_startfix$start.arg, arg_startfix$fix.arg,     argddistname, hasnodefaultval): 'start' must specify names which are arguments to 'distr'.

# (3) fit a discrete distribution (Poisson)
#

set.seed(1234)
x2 <- rpois(n=30,lambda = 2)
qmedist(x2, "pois", probs=1/2)
#> $estimate
#> lambda 
#>    1.7 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.25
#> 
#> $hessian
#>        lambda
#> lambda      0
#> 
#> $optim.function
#> [1] "optim"
#> 
#> $optim.method
#> [1] "BFGS"
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
#>        1        1 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -46.18434
#> 
#> $probs
#> [1] 0.5
#> 

# (4) fit a finite-support distribution (beta)
#

set.seed(1234)
x3 <- rbeta(n=100,shape1=5, shape2=10)
qmedist(x3, "beta", probs=c(1/3, 2/3))
#> $estimate
#>    shape1    shape2 
#>  5.820826 14.053655 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 3.889731e-12
#> 
#> $hessian
#>              shape1        shape2
#> shape1  0.002714767 -0.0010963293
#> shape2 -0.001096329  0.0004477195
#> 
#> $optim.function
#> [1] "optim"
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
#>       89       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] 76.02016
#> 
#> $probs
#> [1] 0.3333333 0.6666667
#> 

# (5) fit frequency distributions on USArrests dataset.
#

x4 <- USArrests$Assault
qmedist(x4, "pois", probs=1/2)
#> $estimate
#> lambda 
#> 170.76 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 144
#> 
#> $hessian
#>        lambda
#> lambda      0
#> 
#> $optim.function
#> [1] "optim"
#> 
#> $optim.method
#> [1] "BFGS"
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
#>        1        1 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1211.705
#> 
#> $probs
#> [1] 0.5
#> 
qmedist(x4, "nbinom", probs=c(1/3, 2/3))
#> $estimate
#>       size         mu 
#>   2.518966 182.313344 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.1111111
#> 
#> $hessian
#>      size mu
#> size    0  0
#> mu      0  0
#> 
#> $optim.function
#> [1] "optim"
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
#>       37       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -292.5969
#> 
#> $probs
#> [1] 0.3333333 0.6666667
#> 
```
