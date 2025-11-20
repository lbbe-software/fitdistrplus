# Fictive survival dataset of a french Male population

100 male individuals randomly taken from `frefictivetable` in
`CASdatasets` package

## Usage

``` r
data(fremale)
```

## Format

`fremale` is a data frame with 3 columns names `AgeIn`, `AgeOut`
respectively for entry age and exit age; `Death` a binary dummy: 1
indicating the death of the individual; 0 a censored observation.

## References

See full dataset `frefictivetable` of `CASdatasets` at
<http://dutangc.perso.math.cnrs.fr/RRepository/>

## Examples

``` r
# (1) load of data
#
data(fremale)
summary(fremale)
#>      AgeIn           AgeOut          Death    
#>  Min.   :23.87   Min.   :30.20   Min.   :0.0  
#>  1st Qu.:47.29   1st Qu.:53.82   1st Qu.:1.0  
#>  Median :63.95   Median :69.49   Median :1.0  
#>  Mean   :60.34   Mean   :67.00   Mean   :0.8  
#>  3rd Qu.:72.00   3rd Qu.:80.23   3rd Qu.:1.0  
#>  Max.   :89.17   Max.   :97.11   Max.   :1.0  
```
