# Pre-fitting procedure

Search good starting values

## Usage

``` r
prefit(data, distr, method = c("mle", "mme", "qme", "mge"), 
  feasible.par, memp=NULL, order=NULL,
  probs=NULL, qtype=7, gof=NULL, fix.arg=NULL, lower, 
  upper, weights=NULL, silent=TRUE, ...)
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
  `"qme"` for 'quantile matching estimation' and `"mge"` for 'maximum
  goodness-of-fit estimation'.

- feasible.par:

  A named list giving the initial values of parameters of the named
  distribution or a function of data computing initial values and
  returning a named list. This argument may be omitted (default) for
  some distributions for which reasonable starting values are computed
  (see the 'details' section of
  [`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md)).
  It may not be into account for closed-form formulas.

- order:

  A numeric vector for the moment order(s). The length of this vector
  must be equal to the number of parameters to estimate.

- memp:

  A function implementing empirical moments, raw or centered but has to
  be consistent with `distr` argument (and `weights` argument).

- probs:

  A numeric vector of the probabilities for which the quantile matching
  is done. The length of this vector must be equal to the number of
  parameters to estimate.

- qtype:

  The quantile type used by the R
  [`quantile`](https://lbbe-software.github.io/fitdistrplus/reference/quantile.md)
  function to compute the empirical quantiles, (default 7 corresponds to
  the default quantile method in R).

- gof:

  A character string coding for the name of the goodness-of-fit distance
  used : "CvM" for Cramer-von Mises distance,"KS" for Kolmogorov-Smirnov
  distance, "AD" for Anderson-Darling distance, "ADR", "ADL", "AD2R",
  "AD2L" and "AD2" for variants of Anderson-Darling distance described
  by Luceno (2006).

- fix.arg:

  An optional named list giving the values of fixed parameters of the
  named distribution or a function of data computing (fixed) parameter
  values and returning a named list. Parameters with fixed value are
  thus NOT estimated by this maximum likelihood procedure. The use of
  this argument is not possible if `method="mme"` and a closed-form
  formula is used.

- weights:

  an optional vector of weights to be used in the fitting process.
  Should be `NULL` or a numeric vector. If non-`NULL`, weighted MLE is
  used, otherwise ordinary MLE.

- silent:

  A logical to remove or show warnings.

- lower:

  Lower bounds on the parameters.

- upper:

  Upper bounds on the parameters.

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

Searching good starting values is achieved by transforming the
parameters (from their constraint interval to the real line) of the
probability distribution. Indeed,

- positive parameters in \\(0,Inf)\\ are transformed using the logarithm
  (typically the scale parameter `sd` of a normal distribution, see
  [Normal](https://rdrr.io/r/stats/Normal.html)),

- parameters in \\(1,Inf)\\ are transformed using the function
  \\log(x-1)\\,

- probability parameters in \\(0,1)\\ are transformed using the logit
  function \\log(x/(1-x))\\ (typically the parameter `prob` of a
  geometric distribution, see
  [Geometric](https://rdrr.io/r/stats/Geometric.html)),

- negative probability parameters in \\(-1,0)\\ are transformed using
  the function \\log(-x/(1+x))\\,

- real parameters are of course not transformed at all, typically the
  `mean` of a normal distribution, see
  [Normal](https://rdrr.io/r/stats/Normal.html).

Once parameters are transformed, an optimization is carried out by a
quasi-Newton algorithm (typically BFGS) and then we transform them back
to original parameter value.

## Value

A named list.

## See also

See
[`mledist`](https://lbbe-software.github.io/fitdistrplus/reference/mledist.md),
[`mmedist`](https://lbbe-software.github.io/fitdistrplus/reference/mmedist.md),
[`qmedist`](https://lbbe-software.github.io/fitdistrplus/reference/qmedist.md),
[`mgedist`](https://lbbe-software.github.io/fitdistrplus/reference/mgedist.md)
for details on parameter estimation. See
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md)
for the main procedure.

## References

Delignette-Muller ML and Dutang C (2015), *fitdistrplus: An R Package
for Fitting Distributions*. Journal of Statistical Software, 64(4),
1-34, [doi:10.18637/jss.v064.i04](https://doi.org/10.18637/jss.v064.i04)
.

## Author

Christophe Dutang and Marie-Laure Delignette-Muller.

## Examples

``` r
set.seed(123) # here just to make random sampling reproducible

# (1) fit of a gamma distribution by maximum likelihood estimation
#
x <- rgamma(1e3, 5/2, 7/2)

prefit(x, "gamma", "mle", list(shape=3, scale=3), lower=-Inf, upper=Inf)
#> $shape
#> [1] 2.671733
#> 
#> $scale
#> [1] 3.901235
#> 
```
