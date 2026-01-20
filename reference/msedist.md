# Maximum spacing estimation of univariate distributions

Fit of univariate distribution by maximizing (log) spacings for non
censored data.

## Usage

``` r
msedist(data, distr, phidiv="KL", power.phidiv=NULL, start = NULL, fix.arg = NULL, 
  optim.method = "default", lower = -Inf, upper = Inf, custom.optim = NULL, 
  weights=NULL, silent = TRUE, gradient = NULL, checkstartfix=FALSE, calcvcov=FALSE, ...)
```

## Arguments

- data:

  A numeric vector for non censored data.

- distr:

  A character string `"name"` naming a distribution for which the
  corresponding quantile function `qname` and the corresponding density
  distribution `dname` must be classically defined.

- phidiv:

  A character string coding for the name of the phi-divergence used :
  `"KL"` for Kullback-Leibler information (corresponds to classic
  maximum spacing estimation), `"J"` for Jeffreys' divergence, `"R"` for
  Renyi's divergence, `"H"` for Hellinger distance, `"V"` for Vajda's
  measure of information, see details.

- power.phidiv:

  If relevant, a numeric for the power used in some phi-divergence :
  should be `NULL` when `phidiv="KL"` or `phidiv="J"` , should be
  positive and different from 1 when `phidiv="R"`, should be greater or
  equal to 1 when `phidiv="H"` or `phidiv="V"`, see details.

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

- weights:

  an optional vector of weights to be used in the fitting process.
  Should be `NULL` or a numeric vector with strictly positive integers
  (typically the number of occurences of each observation). If
  non-`NULL`, weighted MSE is used, otherwise ordinary MSE.

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

The `msedist` function numerically maximizes a phi-divergence function
of spacings, where spacings are the differences of the cumulative
distribution function evaluated at the sorted dataset. The classical
maximum spacing estimation (MSE) was introduced by Cheng and Amin (1986)
and Ranneby (1984) independently where the phi-diverence is the
logarithm, see Anatolyev and Kosenok (2005) for a link between MSE and
maximum likelihood estimation.

MSE was generalized by Ranneby and Ekstrom (1997) by allowing different
phi-divergence function. Generalized MSE maximizes \$\$
S_n(\theta)=\frac{1}{n+1}\sum\_{i=1}^{n+1} \phi\left(F(x\_{(i)};
\theta)-F(x\_{(i-1)}; \theta) \right), \$\$ where \\F(;\theta)\\ is the
parametric distribution function to be fitted, \\\phi\\ is the
phi-divergence function, \\x\_{(1)}\<\dots\<x\_{(n)}\\ is the sorted
sample, \\x\_{(0)}=-\infty\\ and \\x\_{(n+1)}=+\infty\\. The possible
phi-divergence function is

- Kullback-Leibler information (when `phidiv="KL"` and corresponds to
  classical MSE) \$\$\phi(x)=\log(x)\$\$

- Jeffreys' divergence (when `phidiv="J"`) \$\$\phi(x)=(1-x)\log(x)\$\$

- Renyi's divergence (when `phidiv="R"` and `power.phidiv=alpha`)
  \$\$\phi(x)=x^\alpha\times\textrm{sign}(1-\alpha) \textrm{ with }
  \alpha\>0, \alpha\neq 1 \$\$

- Hellinger distance (when `phidiv="H"` and `power.phidiv=p`)
  \$\$\phi(x)=-\|1-x^{1/p}\|^p \textrm{ with } p\ge 1 \$\$

- Vajda's measure of information (when `phidiv="V"` and
  `power.phidiv=beta`) \$\$\phi(x)=-\|1-x\|^\beta \textrm{ with }
  \beta\ge 1 \$\$

The optimization process is the same as
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
see the 'details' section of that function.

This function is not intended to be called directly but is internally
called in
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
and
[`bootdist`](https://lbbe-software.github.io/fitdistrplus/reference/bootdist.md).

This function is intended to be used only with non-censored data.

NB: if your data values are particularly small or large, a scaling may
be needed before the optimization process, see
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)'s
examples.

## Value

`msedist` returns a list with following components,

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

- phidiv:

  The character string coding for the name of the phi-divergence used
  either `"KL"`, `"J"`, `"R"`, `"H"` or `"V"`.

- power.phidiv:

  Either `NULL` or a numeric for the power used in the phi-divergence.

## See also

See
[`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md),
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
[`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md),
[`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md),
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
for other estimation methods.

Please visit the [Frequently Asked
Questions](https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html).

## References

Anatolyev, S., and Kosenok, G. (2005). *An alternative to maximum
likelihood based on spacings*. Econometric Theory, 21(2), 472-476,
[doi:10.1017/S0266466605050255](https://doi.org/10.1017/S0266466605050255)
.

Cheng, R.C.H. and N.A.K. Amin (1983) *Estimating parameters in
continuous univariate distributions with a shifted origin*. Journal of
the Royal Statistical Society Series B 45, 394-403,
[doi:10.1111/j.2517-6161.1983.tb01268.x](https://doi.org/10.1111/j.2517-6161.1983.tb01268.x)
.

Ranneby, B. (1984) *The maximum spacing method: An estimation method
related to the maximum likelihood method*. Scandinavian Journal of
Statistics 11, 93-112.

Ranneby, B. and Ekstroem, M. (1997). *Maximum spacing estimates based on
different metrics*. Umea universitet.

## Author

Marie-Laure Delignette-Muller and Christophe Dutang.

## Examples

``` r
set.seed(123) # here just to make random sampling reproducible

# (1) Fit of a Weibull distribution to serving size data by maximum 
# spacing estimation
# 

data(groundbeef)
serving <- groundbeef$serving
msedist(serving, "weibull")
#> $estimate
#>     shape     scale 
#>  1.423799 80.894950 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 3.789824
#> 
#> $hessian
#>              shape         scale
#> shape  0.792656647 -0.0043440632
#> scale -0.004344063  0.0002995895
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
#>       59       NA 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1287.97
#> 
#> $phidiv
#> [1] "KL"
#> 
#> $power.phidiv
#> NULL
#> 



# (2) Fit of an exponential distribution 
# 

x1 <- rexp(1e3)
#the convergence is quick
msedist(x1, "exp", control=list(trace=0, REPORT=1))
#> $estimate
#>     rate 
#> 0.967625 
#> 
#> $convergence
#> [1] 0
#> 
#> $value
#> [1] 7.516802
#> 
#> $hessian
#>          rate
#> rate 1.066843
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
#>       11        2 
#> 
#> $optim.message
#> NULL
#> 
#> $loglik
#> [1] -1029.544
#> 
#> $phidiv
#> [1] "KL"
#> 
#> $power.phidiv
#> NULL
#> 
```
