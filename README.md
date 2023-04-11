# Help to Fit of a Parametric Distribution to Non-Censored or Censored Data 

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/fitdistrplus)](https://cran.r-project.org/package=fitdistrplus)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/fitdistrplus)](https://cran.r-project.org/package=fitdistrplus)
[![R-CMD-check](https://github.com/aursiber/fitdistrplus/workflows/R-CMD-check/badge.svg)](https://github.com/aursiber/fitdistrplus/actions)

`fitdistrplus` extends the `fitdistr()` function (of the `MASS` package) with several functions to help the fit of a parametric distribution to non-censored or censored data. Censored data may contain left censored, right censored and interval censored values, with several lower and upper bounds. In addition to maximum likelihood estimation (MLE), the package provides moment matching (MME), quantile matching (QME) and maximum goodness-of-fit estimation (MGE) methods (available only for non-censored data). Weighted versions of MLE, MME and QME are available.


## The package

The stable version of `fitdistrplus` can be installed from CRAN using:
```r
install.packages("fitdistrplus")
```

The development version of `fitdistrplus` can be installed from GitHub (`remotes` needed):
```r
if (!requireNamespace("remotes", quietly = TRUE))
   install.packages("remotes")
   
remotes::install_github("aursiber/fitdistrplus")
``` 

Finally load the package in your current R session with the following R command:
```r
library(fitdistrplus)
```

## Documentation

Three **vignettes** are attached to the `fitdistrplus` package:

- <a href="https://aursiber.github.io/fitdistrplus/articles/fitdistrplus_vignette.html" target="_blank">Overview of the fitdistrplus package</a>
- <a href="https://aursiber.github.io/fitdistrplus/articles/Optimalgo.html" target="_blank">Which optimization algorithm to choose?</a>
- <a href="https://aursiber.github.io/fitdistrplus/articles/FAQ.html" target="_blank">Frequently Asked Questions</a>


## Authors & Contacts

Issues can be reported on https://github.com/aursiber/fitdistrplus/issues.

- Marie-Laure Delignette-Muller: marielaure.delignettemuller@vetagro-sup.fr
- Christophe Dutang: dutangc@gmail.com
- Aurélie Siberchicot: aurelie.siberchicot@univ-lyon1.fr


## Citation

If you use `fitdistrplus`, you should cite: <br />
Marie Laure Delignette-Muller, Christophe Dutang (2015). 
*fitdistrplus: An R Package for Fitting Distributions.*
Journal of Statistical Software.
<a href="https://www.jstatsoft.org/article/view/v064i04" target="_blank">https://www.jstatsoft.org/article/view/v064i04</a>
DOI 10.18637/jss.v064.i04.
