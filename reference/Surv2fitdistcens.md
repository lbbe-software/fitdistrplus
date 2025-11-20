# Handling of data formated as in the survival package for use in fitdistcens()

Provide a function to prepare a data frame needed by fitdistcens() from
data classically coded when using the Surv() function of the survival
package

## Usage

``` r
Surv2fitdistcens(time, time2, event,
                      type = c('right', 'left', 'interval', 'interval2'))
```

## Arguments

- time:

  for right censored data, this is the follow up time. For interval
  data, the first argument is the starting time for the interval.

- event:

  The status indicator, normally `0`=alive, `1`=dead. Other choices are
  `TRUE/FALSE` (`TRUE` = death) or `1/2` (`2`=death). For interval
  censored data, the status indicator is `0`=right censored, `1`=event
  at time, `2`=left censored, `3`=interval censored. For factor data,
  assume that it has only two levels with the second level coding death.

- time2:

  ending time of the interval for interval censored. Intervals are
  assumed to be open on the left and closed on the right, (start, end\].

- type:

  character string specifying the type of censoring. Possible values are
  `"right"`, `"left"`, `"interval"`, `"interval2"`.

## Details

`Surv2fitdistcens` makes a `data.frame` with two columns respectively
named `left` and `right`, describing each observed value as an interval
as required in fitdistcens(): the `left` column contains either `NA` for
left-censored observations, the left bound of the interval for
interval-censored observations, or the observed value for non-censored
observations. The right column contains either `NA` for right-censored
observations, the right bound of the interval for interval censored
observations, or the observed value for non-censored observations.

## Value

`Surv2fitdistcens` returns a data.frame with two columns respectively
named `left` and `right`.

## See also

See
[`fitdistrplus`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistrplus.md)
for an overview of the package. See
[`fitdistcens`](https://lbbe-software.github.io/fitdistrplus/reference/fitdistcens.md)
for fitting of univariate distributions to censored data and
[`fremale`](https://lbbe-software.github.io/fitdistrplus/reference/fremale.md)
for the full dataset used in examples below. See
[`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) for survival
objects which use the same arguments.

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
# (1) randomized fictive survival data - right-censored
#
origdata <- data.frame(rbind(
c(   43.01, 55.00,     0),
c(   36.37, 47.17,     0),
c(   33.10, 34.51,     0),
c(   71.00, 81.15,     1),
c(   80.89, 81.91,     1),
c(   67.81, 78.48,     1),
c(   73.98, 76.92,     1),
c(   53.19, 54.80,     1)))
colnames(origdata) <- c("AgeIn", "AgeOut", "Death")

# add of follow-up time (for type = "right" in Surv())
origdata$followuptime <- origdata$AgeOut - origdata$AgeIn
origdata
#>   AgeIn AgeOut Death followuptime
#> 1 43.01  55.00     0        11.99
#> 2 36.37  47.17     0        10.80
#> 3 33.10  34.51     0         1.41
#> 4 71.00  81.15     1        10.15
#> 5 80.89  81.91     1         1.02
#> 6 67.81  78.48     1        10.67
#> 7 73.98  76.92     1         2.94
#> 8 53.19  54.80     1         1.61

### use of default survival type "right"
# in Surv()
survival::Surv(time = origdata$followuptime, event = origdata$Death, type = "right")
#> [1] 11.99+ 10.80+  1.41+ 10.15   1.02  10.67   2.94   1.61 
# for fitdistcens()
Surv2fitdistcens(origdata$followuptime, event = origdata$Death, type = "right")
#>    left right
#> 1 11.99    NA
#> 2 10.80    NA
#> 3  1.41    NA
#> 4 10.15 10.15
#> 5  1.02  1.02
#> 6 10.67 10.67
#> 7  2.94  2.94
#> 8  1.61  1.61

# use of survival type "interval" 
# in Surv()
survival::Surv(time = origdata$followuptime, time2 = origdata$followuptime, 
          event = origdata$Death, type = "interval")
#> [1] 11.99+ 10.80+  1.41+ 10.15   1.02  10.67   2.94   1.61 
# for fitdistcens()
Surv2fitdistcens(time = origdata$followuptime, time2 = origdata$followuptime, 
          event = origdata$Death, type = "interval") 
#>    left right
#> 1 11.99    NA
#> 2 10.80    NA
#> 3  1.41    NA
#> 4 10.15 10.15
#> 5  1.02  1.02
#> 6 10.67 10.67
#> 7  2.94  2.94
#> 8  1.61  1.61

# use of survival type "interval2" 
origdata$survivalt1 <- origdata$followuptime
origdata$survivalt2 <- origdata$survivalt1
origdata$survivalt2[1:3] <- Inf
origdata
#>   AgeIn AgeOut Death followuptime survivalt1 survivalt2
#> 1 43.01  55.00     0        11.99      11.99        Inf
#> 2 36.37  47.17     0        10.80      10.80        Inf
#> 3 33.10  34.51     0         1.41       1.41        Inf
#> 4 71.00  81.15     1        10.15      10.15      10.15
#> 5 80.89  81.91     1         1.02       1.02       1.02
#> 6 67.81  78.48     1        10.67      10.67      10.67
#> 7 73.98  76.92     1         2.94       2.94       2.94
#> 8 53.19  54.80     1         1.61       1.61       1.61
survival::Surv(time = origdata$survivalt1, time2 = origdata$survivalt2, 
type = "interval2")
#> [1] 11.99+ 10.80+  1.41+ 10.15   1.02  10.67   2.94   1.61 
Surv2fitdistcens(origdata$survivalt1, time2 = origdata$survivalt2, 
                type = "interval2")
#>    left right
#> 1 11.99    NA
#> 2 10.80    NA
#> 3  1.41    NA
#> 4 10.15 10.15
#> 5  1.02  1.02
#> 6 10.67 10.67
#> 7  2.94  2.94
#> 8  1.61  1.61


# (2) Other examples with various left, right and interval censored values
#
# with left censored data
(d1 <- data.frame(time = c(2, 5, 3, 7), ind = c(0, 1, 1, 1)))
#>   time ind
#> 1    2   0
#> 2    5   1
#> 3    3   1
#> 4    7   1
survival::Surv(time = d1$time, event = d1$ind, type = "left")
#> [1] 2- 5  3  7 
Surv2fitdistcens(time = d1$time, event = d1$ind, type = "left")
#>   left right
#> 1   NA     2
#> 2    5     5
#> 3    3     3
#> 4    7     7

(d1bis <- data.frame(t1 = c(2, 5, 3, 7), t2 = c(2, 5, 3, 7), 
  censtype = c(2, 1, 1, 1)))
#>   t1 t2 censtype
#> 1  2  2        2
#> 2  5  5        1
#> 3  3  3        1
#> 4  7  7        1
survival::Surv(time = d1bis$t1, time2 = d1bis$t2, 
  event = d1bis$censtype, type = "interval")
#> [1] 2- 5  3  7 
Surv2fitdistcens(time = d1bis$t1, time2 = d1bis$t2, 
  event = d1bis$censtype, type = "interval")
#>   left right
#> 1   NA     2
#> 2    5     5
#> 3    3     3
#> 4    7     7

# with interval, left and right censored data
(d2 <- data.frame(t1 = c(-Inf, 2, 3, 4, 3, 7), t2 = c(2, 5, 3, 7, 8, Inf)))
#>     t1  t2
#> 1 -Inf   2
#> 2    2   5
#> 3    3   3
#> 4    4   7
#> 5    3   8
#> 6    7 Inf
survival::Surv(time = d2$t1, time2 = d2$t2, type = "interval2")
#> [1] 2-     [2, 5] 3      [4, 7] [3, 8] 7+    
Surv2fitdistcens(time = d2$t1, time2 = d2$t2, type = "interval2")
#>   left right
#> 1   NA     2
#> 2    2     5
#> 3    3     3
#> 4    4     7
#> 5    3     8
#> 6    7    NA

(d2bis <- data.frame(t1 = c(2, 2, 3, 4, 3, 7), t2 = c(2, 5, 3, 7, 8, 7), 
  censtype = c(2,3,1,3,3,0)))
#>   t1 t2 censtype
#> 1  2  2        2
#> 2  2  5        3
#> 3  3  3        1
#> 4  4  7        3
#> 5  3  8        3
#> 6  7  7        0
survival::Surv(time = d2bis$t1, time2 = d2bis$t2, 
  event = d2bis$censtype, type = "interval")
#> [1] 2-     [2, 5] 3      [4, 7] [3, 8] 7+    
Surv2fitdistcens(time = d2bis$t1, time2 = d2bis$t2, 
  event = d2bis$censtype, type = "interval")
#>   left right
#> 1   NA     2
#> 2    2     5
#> 3    3     3
#> 4    4     7
#> 5    3     8
#> 6    7    NA
```
