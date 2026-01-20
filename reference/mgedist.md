# Maximum goodness-of-fit fit of univariate continuous distributions

Fit of univariate continuous distribution by maximizing goodness-of-fit
(or minimizing distance) for non censored data.

## Usage

``` r
mgedist(data, distr, gof = "CvM", start = NULL, fix.arg = NULL, optim.method = "default", 
  lower = -Inf, upper = Inf, custom.optim = NULL, silent = TRUE, gradient = NULL, 
  checkstartfix=FALSE, calcvcov=FALSE, ...)
```

## Arguments

- data:

  A numeric vector for non censored data.

- distr:

  A character string `"name"` naming a distribution for which the
  corresponding quantile function `qname` and the corresponding density
  distribution `dname` must be classically defined.

- gof:

  A character string coding for the name of the goodness-of-fit distance
  used : `"CvM"` for Cramer-von Mises distance, `"KS"` for
  Kolmogorov-Smirnov distance, `"AD"` for Anderson-Darling distance,
  `"ADR"`, `"ADL"`, `"AD2R"`, `"AD2L"` and `"AD2"` for variants of
  Anderson-Darling distance described by Luceno (2006).

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

  a function carrying the optimization.

- silent:

  A logical to remove or show warnings when bootstraping.

- gradient:

  A function to return the gradient of the gof distance for the
  `"BFGS"`, `"CG"` and `"L-BFGS-B"` methods. If it is `NULL`, a
  finite-difference approximation will be used.

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

The `mgedist` function numerically maximizes goodness-of-fit, or
minimizes a goodness-of-fit distance coded by the argument `gof`. One
may use one of the classical distances defined in Stephens (1986), the
Cramer-von Mises distance (`"CvM"`), the Kolmogorov-Smirnov distance
(`"KS"`) or the Anderson-Darling distance (`"AD"`) which gives more
weight to the tails of the distribution, or one of the variants of this
last distance proposed by Luceno (2006). The right-tail AD (`"ADR"`)
gives more weight only to the right tail, the left-tail AD (`"ADL"`)
gives more weight only to the left tail. Either of the tails, or both of
them, can receive even larger weights by using second order
Anderson-Darling Statistics (using `"AD2R"`, `"AD2L"` or `"AD2"`).

The optimization process is the same as
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
see the 'details' section of that function.

This function is not intended to be called directly but is internally
called in
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
and
[`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md).

This function is intended to be used only with continuous distributions
and weighted maximum goodness-of-fit estimation is not allowed.

NB: if your data values are particularly small or large, a scaling may
be needed before the optimization process. See example (4).

## Value

`mgedist` returns a list with following components,

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

- gof:

  the code of the goodness-of-fit distance maximized.

## See also

See
[`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md),
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
[`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md),
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
for other estimation methods.

Please visit the [Frequently Asked
Questions](https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html).

## References

Luceno A (2006), *Fitting the generalized Pareto distribution to data
using maximum goodness-of-fit estimators*. Computational Statistics and
Data Analysis, 51, 904-917,
[doi:10.1016/j.csda.2005.09.011](https://doi.org/10.1016/j.csda.2005.09.011)
.

Stephens MA (1986), *Tests based on edf statistics*. In Goodness-of-fit
techniques (D'Agostino RB and Stephens MA, eds), Marcel Dekker, New
York, pp. 97-194.

Delignette-Muller ML and Dutang C (2015), *fitdistrplus: An R Package
for Fitting Distributions*. Journal of Statistical Software, 64(4),
1-34, [doi:10.18637/jss.v064.i04](https://doi.org/10.18637/jss.v064.i04)
.

## Author

Marie-Laure Delignette-Muller and Christophe Dutang.

## Examples

``` r
set.seed(123) # here just to make random sampling reproducible

# (1) Fit of a Weibull distribution to serving size data by maximum 
# goodness-of-fit estimation using all the distances available
# 

data(groundbeef)
serving <- groundbeef$serving
mgedist(serving, "weibull", gof="CvM")
#> $estimate
#>     shape     scale 
#>  2.093204 82.660014 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.002581367
#> 
#> $hessian
#>              shape        scale
#> shape 0.0159565105 3.639558e-04
#> scale 0.0003639558 9.522745e-05
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
#>       65       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1255.623
#> 
#> $gof
#> [1] "CvM"
#> 
mgedist(serving, "weibull", gof="KS")
#> $estimate
#>     shape     scale 
#>  2.065634 81.450487 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.112861
#> 
#> $hessian
#>            shape    scale
#> shape 122.668263 6.509057
#> scale   6.509057 7.599584
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
#>      127       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1255.975
#> 
#> $gof
#> [1] "KS"
#> 
mgedist(serving, "weibull", gof="AD")
#> $estimate
#>     shape     scale 
#>  2.125425 82.890502 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.0137836
#> 
#> $hessian
#>              shape        scale
#> shape 0.1158157367 0.0007180241
#> scale 0.0007180241 0.0005332051
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
#>       67       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1255.393
#> 
#> $gof
#> [1] "AD"
#> 
mgedist(serving, "weibull", gof="ADR")
#> $estimate
#>     shape     scale 
#>  2.072035 82.762593 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.006340469
#> 
#> $hessian
#>              shape         scale
#> shape  0.053243854 -0.0013083937
#> scale -0.001308394  0.0003140377
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
#>       69       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1255.837
#> 
#> $gof
#> [1] "ADR"
#> 
mgedist(serving, "weibull", gof="ADL")
#> $estimate
#>     shape     scale 
#>  2.197498 82.016005 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.007267475
#> 
#> $hessian
#>             shape        scale
#> shape 0.060343316 0.0021420124
#> scale 0.002142012 0.0002184993
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
#>       65       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1255.415
#> 
#> $gof
#> [1] "ADL"
#> 
mgedist(serving, "weibull", gof="AD2R")
#> $estimate
#>    shape    scale 
#>  1.90328 81.33464 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.04552816
#> 
#> $hessian
#>             shape        scale
#> shape  1.31736538 -0.041034447
#> scale -0.04103445  0.002056365
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
#>       69       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1259.112
#> 
#> $gof
#> [1] "AD2R"
#> 
mgedist(serving, "weibull", gof="AD2L")
#> $estimate
#>     shape     scale 
#>  2.483836 78.252113 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.0385314
#> 
#> $hessian
#>            shape        scale
#> shape 0.44689737 0.0161843919
#> scale 0.01618439 0.0009217762
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
#>       69       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1265.933
#> 
#> $gof
#> [1] "AD2L"
#> 
mgedist(serving, "weibull", gof="AD2")
#> $estimate
#>     shape     scale 
#>  2.081168 85.281194 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.1061089
#> 
#> $hessian
#>             shape       scale
#> shape  2.10614403 -0.04170905
#> scale -0.04170905  0.00299467
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
#>       69       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1256.313
#> 
#> $gof
#> [1] "AD2"
#> 


# (2) Fit of a uniform distribution using Cramer-von Mises or
# Kolmogorov-Smirnov distance
# 

u <- runif(100,min=5,max=10)
mgedist(u,"unif",gof="CvM")
#> $estimate
#>      min      max 
#> 5.037213 9.948626 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.0002131867
#> 
#> $hessian
#>            min        max
#> min 0.02679927 0.01382797
#> max 0.01382797 0.02679853
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
#>       79       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -Inf
#> 
#> $gof
#> [1] "CvM"
#> 
mgedist(u,"unif",gof="KS")
#> $estimate
#>      min      max 
#> 5.015314 9.888292 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.0428921
#> 
#> $hessian
#>           min       max
#> min 145.68231  19.49759
#> max  19.49759 142.63348
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
#>      113       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -Inf
#> 
#> $gof
#> [1] "KS"
#> 

# (3) Fit of a triangular distribution using Cramer-von Mises or
# Kolmogorov-Smirnov distance
# 

# \donttest{
require("mc2d")
t <- rtriang(100,min=5,mode=6,max=10)
mgedist(t,"triang",start = list(min=4, mode=6,max=9),gof="CvM")
#> Warning: Some parameter names have no starting/fixed value but have a default value: mean.
#> $estimate
#>      min     mode      max 
#> 5.429906 5.971934 9.636604 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.0003794627
#> 
#> $hessian
#>             min       mode        max
#> min  0.03331760 0.03504449 0.01591392
#> mode 0.03504449 0.03908299 0.01802510
#> max  0.01591392 0.01802510 0.01558397
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
#>       84       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -Inf
#> 
#> $gof
#> [1] "CvM"
#> 
mgedist(t,"triang",start = list(min=4, mode=6,max=9),gof="KS")
#> Warning: Some parameter names have no starting/fixed value but have a default value: mean.
#> $estimate
#>      min     mode      max 
#> 5.529346 5.861222 9.676341 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 0.04535746
#> 
#> $hessian
#>            min      mode       max
#> min  138.62266 138.63822  48.79111
#> mode 138.63822 150.67871  50.89331
#> max   48.79111  50.89331 113.65007
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
#>      206       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -Inf
#> 
#> $gof
#> [1] "KS"
#> 
# }

# (4) scaling problem
# the simulated dataset (below) has particularly small values, hence without scaling (10^0),
# the optimization raises an error. The for loop shows how scaling by 10^i
# for i=1,...,6 makes the fitting procedure work correctly.

x2 <- rnorm(100, 1e-4, 2e-4)
for(i in 6:0)
    cat(i, try(mgedist(x*10^i,"cauchy")$estimate, silent=TRUE), "\n")
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
```
