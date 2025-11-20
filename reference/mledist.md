# Maximum likelihood fit of univariate distributions

Fit of univariate distributions using maximum likelihood for censored or
non censored data.

## Usage

``` r
mledist(data, distr, start = NULL, fix.arg = NULL, optim.method = "default", 
    lower = -Inf, upper = Inf, custom.optim = NULL, weights = NULL, silent = TRUE, 
    gradient = NULL, checkstartfix=FALSE, calcvcov=FALSE, ...)
```

## Arguments

- data:

  A numeric vector for non censored data or a dataframe of two columns
  respectively named `left` and `right`, describing each observed value
  as an interval for censored data. In that case the `left` column
  contains either `NA` for left censored observations, the left bound of
  the interval for interval censored observations, or the observed value
  for non-censored observations. The `right` column contains either `NA`
  for right censored observations, the right bound of the interval for
  interval censored observations, or the observed value for non-censored
  observations.

- distr:

  A character string `"name"` naming a distribution for which the
  corresponding density function `dname` and the corresponding
  distribution function `pname` must be classically defined.

- start:

  A named list giving the initial values of parameters of the named
  distribution or a function of data computing initial values and
  returning a named list. This argument may be omitted (default) for
  some distributions for which reasonable starting values are computed
  (see details).

- fix.arg:

  An optional named list giving the values of fixed parameters of the
  named distribution or a function of data computing (fixed) parameter
  values and returning a named list. Parameters with fixed value are
  thus NOT estimated by this maximum likelihood procedure.

- optim.method:

  `"default"` (see details) or an optimization method to pass to
  [`optim`](https://rdrr.io/r/stats/optim.html).

- lower:

  Left bounds on the parameters for the `"L-BFGS-B"` method (see
  [`optim`](https://rdrr.io/r/stats/optim.html)).

- upper:

  Right bounds on the parameters for the `"L-BFGS-B"` method (see
  [`optim`](https://rdrr.io/r/stats/optim.html)).

- custom.optim:

  a function carrying the MLE optimisation (see details).

- weights:

  an optional vector of weights to be used in the fitting process.
  Should be `NULL` or a numeric vector with strictly positive integers
  (typically the number of occurences of each observation). If
  non-`NULL`, weighted MLE is used, otherwise ordinary MLE.

- silent:

  A logical to remove or show warnings when bootstraping.

- gradient:

  A function to return the gradient of the log-likelihood for the
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

This function is not intended to be called directly but is internally
called in
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
and
[`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md)
when used with the maximum likelihood method and
[`fitdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.md)
and
[`bootdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/bootdistcens.md).

It is assumed that the `distr` argument specifies the distribution by
the probability density function and the cumulative distribution
function (d, p). The quantile function and the random generator function
(q, r) may be needed by other function such as
[`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md),
[`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md),
[`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md),
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md),[`fitdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.md),
[`bootdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/bootdistcens.md)
and
[`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md).

For the following named distributions, reasonable starting values will
be computed if `start` is omitted (i.e. `NULL`) : `"norm"`, `"lnorm"`,
`"exp"` and `"pois"`, `"cauchy"`, `"gamma"`, `"logis"`, `"nbinom"`
(parametrized by mu and size), `"geom"`, `"beta"`, `"weibull"` from the
`stats` package; all distributions (except phase-type distributions)
from the `actuar` package. Note that these starting values may not be
good enough if the fit is poor. The function uses a closed-form formula
to fit the uniform distribution. If `start` is a list, then it should be
a named list with the same names as in the d,p,q,r functions of the
chosen distribution. If `start` is a function of data, then the function
should return a named list with the same names as in the d,p,q,r
functions of the chosen distribution.

The `mledist` function allows user to set a fixed values for some
parameters. As for `start`, if `fix.arg` is a list, then it should be a
named list with the same names as in the d,p,q,r functions of the chosen
distribution. If `fix.arg` is a function of data, then the function
should return a named list with the same names as in the d,p,q,r
functions of the chosen distribution.

When `custom.optim=NULL` (the default), maximum likelihood estimations
of the distribution parameters are computed with the R base
[`optim`](https://rdrr.io/r/stats/optim.html) or
[`constrOptim`](https://rdrr.io/r/stats/constrOptim.html). If no finite
bounds (`lower=-Inf` and `upper=Inf`) are supplied,
[`optim`](https://rdrr.io/r/stats/optim.html) is used with the method
specified by `optim.method`. Note that `optim.method="default"` means
`optim.method="Nelder-Mead"` for distributions with at least two
parameters and `optim.method="BFGS"` for distributions with only one
parameter. If finite bounds are supplied (among `lower` and `upper`) and
`gradient != NULL`,
[`constrOptim`](https://rdrr.io/r/stats/constrOptim.html) is used. If
finite bounds are supplied (among `lower` and `upper`) and
`gradient == NULL`,
[`constrOptim`](https://rdrr.io/r/stats/constrOptim.html) is used when
`optim.method="Nelder-Mead"`;
[`optim`](https://rdrr.io/r/stats/optim.html) is used when
`optim.method="L-BFGS-B"` or `"Brent"`; in other case, an error is
raised (same behavior as
[`constrOptim`](https://rdrr.io/r/stats/constrOptim.html)).

When errors are raised by [`optim`](https://rdrr.io/r/stats/optim.html),
it's a good idea to start by adding traces during the optimization
process by adding `control=list(trace=1, REPORT=1)`.

If `custom.optim` is not `NULL`, then the user-supplied function is used
instead of the R base [`optim`](https://rdrr.io/r/stats/optim.html). The
`custom.optim` must have (at least) the following arguments `fn` for the
function to be optimized, `par` for the initialized parameters.
Internally the function to be optimized will also have other arguments,
such as `obs` with observations and `ddistname` with distribution name
for non censored data (Beware of potential conflicts with optional
arguments of `custom.optim`). It is assumed that `custom.optim` should
carry out a MINIMIZATION. Finally, it should return at least the
following components `par` for the estimate, `convergence` for the
convergence code, `value` for `fn(par)`, `hessian`, `counts` for the
number of calls (function and gradient) and `message` (default to
`NULL`) for the error message when `custom.optim` raises an error, see
the returned value of [`optim`](https://rdrr.io/r/stats/optim.html). See
examples in
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
and
[`fitdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.md).

Optionally, a vector of `weights` can be used in the fitting process. By
default (when `weigths=NULL`), ordinary MLE is carried out, otherwise
the specified weights are used to balance the log-likelihood
contributions. It is not yet possible to take into account weights in
functions `plotdist`, `plotdistcens`, `plot.fitdist`,
`plot.fitdistcens`, `cdfcomp`, `cdfcompcens`, `denscomp`, `ppcomp`,
`qqcomp`, `gofstat`, `descdist`, `bootdist`, `bootdistcens` and
`mgedist`. (developments planned in the future).

NB: if your data values are particularly small or large, a scaling may
be needed before the optimization process. See Example (7).

## Value

`mledist` returns a list with following components,

- estimate:

  the parameter estimates.

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

- value:

  the minimal value reached for the criterion to minimize.

- hessian:

  a symmetric matrix computed by
  [`optim`](https://rdrr.io/r/stats/optim.html) as an estimate of the
  Hessian at the solution found or computed in the user-supplied
  optimization function. It is used in `fitdist` to estimate standard
  errors.

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

- method:

  `"closed formula"` if appropriate otherwise `NULL`.

## See also

See
[`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md),
[`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md),
[`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md),
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md),[`fitdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.md)
for other estimation methods,
[`optim`](https://rdrr.io/r/stats/optim.html),
[`constrOptim`](https://rdrr.io/r/stats/constrOptim.html) for
optimization routines,
[`bootdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/bootdistcens.md)
and
[`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md)
for bootstrap, and
[`llplot`](https://lbbe-software.github.io/fitdistrplus/reference/logLik-plot.md)
for plotting the (log)likelihood.

Please visit the [Frequently Asked
Questions](https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html).

## References

Venables WN and Ripley BD (2002), *Modern applied statistics with S*.
Springer, New York, pp. 435-446,
[doi:10.1007/978-0-387-21706-2](https://doi.org/10.1007/978-0-387-21706-2)
.

Delignette-Muller ML and Dutang C (2015), *fitdistrplus: An R Package
for Fitting Distributions*. Journal of Statistical Software, 64(4),
1-34, [doi:10.18637/jss.v064.i04](https://doi.org/10.18637/jss.v064.i04)
.

## Author

Marie-Laure Delignette-Muller and Christophe Dutang.

## Examples

``` r
# (1) basic fit of a normal distribution with maximum likelihood estimation
#

set.seed(1234)
x1 <- rnorm(n=100)
mledist(x1,"norm")
#> $estimate
#>       mean         sd 
#> -0.1567617  0.9993707 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 1.418309
#> 
#> $hessian
#>               mean            sd
#> mean  1.001260e+00 -5.551115e-11
#> sd   -5.551115e-11  2.002538e+00
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
#>       43       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -141.8309
#> 
#> $vcov
#> NULL
#> 

# (2) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view dedicated to probability distributions

dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
mledist(x1,"gumbel",start=list(a=10,b=5))
#> Error in checkparamlist(arg_startfix$start.arg, arg_startfix$fix.arg,     argddistname, hasnodefaultval): 'start' must specify names which are arguments to 'distr'.

# (3) fit of a discrete distribution (Poisson)
#

set.seed(1234)
x2 <- rpois(n=30,lambda = 2)
mledist(x2,"pois")
#> $estimate
#> lambda 
#>    1.7 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 1.539478
#> 
#> $hessian
#>           lambda
#> lambda 0.5882357
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
#>        4        1 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -46.18434
#> 
#> $vcov
#> NULL
#> 

# (4) fit a finite-support distribution (beta)
#

set.seed(1234)
x3 <- rbeta(n=100,shape1=5, shape2=10)
mledist(x3,"beta")
#> $estimate
#>    shape1    shape2 
#>  4.859798 10.918841 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] -0.7833052
#> 
#> $hessian
#>             shape1      shape2
#> shape1  0.16295311 -0.06542752
#> shape2 -0.06542752  0.03047900
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
#>       47       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] 78.33052
#> 
#> $vcov
#> NULL
#> 


# (5) fit frequency distributions on USArrests dataset.
#

x4 <- USArrests$Assault
mledist(x4, "pois")
#> $estimate
#> lambda 
#> 170.76 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 24.2341
#> 
#> $hessian
#>             lambda
#> lambda 0.005856174
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
#>        2        1 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1211.705
#> 
#> $vcov
#> NULL
#> 
mledist(x4, "nbinom")
#> $estimate
#>       size         mu 
#>   3.822579 170.747853 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 5.806593
#> 
#> $hessian
#>               size            mu
#> size  3.518616e-02 -3.985701e-07
#> mu   -3.985701e-07  1.282598e-04
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
#>       47       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -290.3297
#> 
#> $vcov
#> NULL
#> 

# (6) fit a continuous distribution (Gumbel) to censored data.
#

data(fluazinam)
log10EC50 <-log10(fluazinam)
# definition of the Gumbel distribution
dgumbel  <-  function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel  <-  function(q,a,b) exp(-exp((a-q)/b))
qgumbel  <-  function(p,a,b) a-b*log(-log(p))

mledist(log10EC50,"gumbel",start=list(a=0,b=2),optim.method="Nelder-Mead")
#> Error in checkparamlist(arg_startfix$start.arg, arg_startfix$fix.arg,     argddistname, hasnodefaultval): 'start' must specify names which are arguments to 'distr'.

# (7) scaling problem
# the simulated dataset (below) has particularly small values, 
# hence without scaling (10^0),
# the optimization raises an error. The for loop shows how scaling by 10^i
# for i=1,...,6 makes the fitting procedure work correctly.

set.seed(1234)
x2 <- rnorm(100, 1e-4, 2e-4)
for(i in 6:0)
    cat(i, try(mledist(x*10^i, "cauchy")$estimate, silent=TRUE), "\n")
#> 6 Error in eval(expr, envir) : object 'x' not found
#>  
#> 5 Error in eval(expr, envir) : object 'x' not found
#>  
#> 4 Error in eval(expr, envir) : object 'x' not found
#>  
#> 3 Error in eval(expr, envir) : object 'x' not found
#>  
#> 2 Error in eval(expr, envir) : object 'x' not found
#>  
#> 1 Error in eval(expr, envir) : object 'x' not found
#>  
#> 0 Error in eval(expr, envir) : object 'x' not found
#>  
        
 
# (17) small example for the zero-modified geometric distribution
#

dzmgeom <- function(x, p1, p2) p1 * (x == 0) + (1-p1)*dgeom(x-1, p2) #pdf
x2 <- c(2,  4,  0, 40,  4, 21,  0,  0,  0,  2,  5,  0,  0, 13,  2) #simulated dataset
initp1 <- function(x) list(p1=mean(x == 0)) #init as MLE
mledist(x2, "zmgeom", fix.arg=initp1, start=list(p2=1/2))
#> Error in checkparamlist(arg_startfix$start.arg, arg_startfix$fix.arg,     argddistname, hasnodefaultval): 'start' must specify names which are arguments to 'distr'.
```
