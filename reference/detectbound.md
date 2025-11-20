# Detect bounds for density function

Manual detection of bounds of parameter of a density function/

## Usage

``` r
detectbound(distname, vstart, obs, fix.arg=NULL, echo=FALSE)
```

## Arguments

- distname:

  A character string `"name"` naming a distribution for which the
  corresponding density function `dname` must be classically defined.

- vstart:

  A named vector giving the initial values of parameters of the named
  distribution.

- obs:

  A numeric vector for non censored data.

- fix.arg:

  An optional named vector giving the values of fixed parameters of the
  named distribution. Default to `NULL`.

- echo:

  A logical to show some traces.

## Details

This function manually tests the following bounds : -1, 0, and 1.

## Value

`detectbound` returns a 2-row matrix with the lower bounds in the first
row and the upper bounds in the second row.

## See also

See
[`fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.md).

Please visit the [Frequently Asked
Questions](https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html).

## References

Delignette-Muller ML and Dutang C (2015), *fitdistrplus: An R Package
for Fitting Distributions*. Journal of Statistical Software, 64(4),
1-34, [doi:10.18637/jss.v064.i04](https://doi.org/10.18637/jss.v064.i04)
.

## Author

Christophe Dutang and Marie-Laure Delignette-Muller.

## Examples

``` r
# case where the density returns a Not-an-Numeric value.
detectbound("exp", c(rate=3), 1:10)
#>      rate
#> lowb    0
#> uppb  Inf
detectbound("binom", c(size=3, prob=1/2), 1:10)
#>      size prob
#> lowb -Inf    0
#> uppb  Inf    1
detectbound("nbinom", c(size=3, prob=1/2), 1:10)
#>      size prob   mu
#> lowb    0    0 -Inf
#> uppb  Inf    1  Inf
```
