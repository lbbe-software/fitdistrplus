# Parasite abundance in insular feral cats

Toxocara cati abundance in feral cats living on Kerguelen island.

## Usage

``` r
data(toxocara)
```

## Format

`toxocara` is a data frame with 1 column (number: number of parasites in
digestive tract)

## Source

Fromont, E., Morvilliers, L., Artois, M., Pontier, D. 2001. Parasite
richness and abundance in insular and mainland feral cats.
*Parasitology*, **123**, 143-151.

## Examples

``` r
# (1) load of data
#
data(toxocara)

# (2) description and plot of data
#
number <- toxocara$number
descdist(number, discrete = TRUE, boot = 11)

#> summary statistics
#> ------
#> min:  0   max:  75 
#> median:  2 
#> mean:  8.679245 
#> estimated sd:  14.29332 
#> estimated skewness:  2.630609 
#> estimated kurtosis:  11.4078 
plotdist(number, discrete = TRUE)


# (3) fit of a Poisson distribution to data
#
fitp <- fitdist(number, "pois")
summary(fitp)
#> Fitting of the distribution ' pois ' by maximum likelihood 
#> Parameters : 
#>        estimate Std. Error
#> lambda 8.679245  0.4046719
#> Loglikelihood:  -507.5334   AIC:  1017.067   BIC:  1019.037 
plot(fitp)


# (4) fit of a negative binomial distribution to data
#
fitnb <- fitdist(number, "nbinom")
summary(fitnb)
#> Fitting of the distribution ' nbinom ' by maximum likelihood 
#> Parameters : 
#>       estimate Std. Error
#> size 0.3971457 0.08289027
#> mu   8.6802520 1.93501002
#> Loglikelihood:  -159.3441   AIC:  322.6882   BIC:  326.6288 
#> Correlation matrix:
#>              size           mu
#> size  1.000000000 -0.000103854
#> mu   -0.000103854  1.000000000
#> 
plot(fitnb)
```
