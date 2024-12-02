# Help to Fit of a Parametric Distribution to Non-Censored or Censored Data  <img src="man/figures/fitdistrplus_hex.png" align="right" style="float:right; height:200px;"/>

[![CRAN_Release_Badge](https://www.r-pkg.org/badges/version-ago/fitdistrplus)](https://cran.r-project.org/package=fitdistrplus)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/fitdistrplus)](https://cran.r-project.org/package=fitdistrplus)
[![R-CMD-check](https://github.com/lbbe-software/fitdistrplus/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lbbe-software/fitdistrplus/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---------------------------

**Please note! Since January 2024, this repository has belonged to the *lbbe-software* organization.**
To avoid confusion, we strongly recommend updating any existing local clones to point to the new 
repository URL. You can do this by using `git remote` on the command line:

`git remote set-url origin git@github.com:lbbe-software/fitdistrplus.git`

or 

`git remote set-url origin https://github.com/lbbe-software/fitdistrplus.git`

---------------------------

`fitdistrplus` extends the `fitdistr()` function (of the `MASS` package) with several functions to help the fit of a parametric distribution to non-censored or censored data. Censored data may contain left censored, right censored and interval censored values, with several lower and upper bounds. In addition to maximum likelihood estimation (MLE), the package provides moment matching (MME), quantile matching (QME) and maximum goodness-of-fit estimation (MGE) methods (available only for non-censored data). Weighted versions of MLE, MME and QME are available.

`fitdistrplus` allows to fit any probability distribution provided by the user and not restricted to base R distributions (see `?Distributions`). We strongly encourage users to visit the CRAN task view 
on [Distributions](https://cran.r-project.org/view=Distributions) proposed by Dutang, Kiener & Swihart (2024).

## The package

The stable version of `fitdistrplus` can be installed from CRAN using:
```r
install.packages("fitdistrplus")
```

The development version of `fitdistrplus` can be installed from GitHub (`remotes` needed):
```r
if (!requireNamespace("remotes", quietly = TRUE))
   install.packages("remotes")
   
remotes::install_github("lbbe-software/fitdistrplus")
``` 

Finally load the package in your current R session with the following R command:
```r
require("fitdistrplus")
```

## Documentation

Four **vignettes** are attached to the `fitdistrplus` package. Two of them are for beginners

- <a href="https://lbbe-software.github.io/fitdistrplus/articles/fitdistrplus_vignette.html" target="_blank">Overview of the fitdistrplus package</a>
- <a href="https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html" target="_blank">Frequently Asked Questions</a>

The last two vignettes deal with advanced topics

- <a href="https://lbbe-software.github.io/fitdistrplus/articles/Optimalgo.html" target="_blank">Which optimization algorithm to choose?</a>
- <a href="https://lbbe-software.github.io/fitdistrplus/articles/starting-values.html" target="_blank">Starting values used in fitdistrplus</a>



## Authors & Contacts

Please read the <a href="https://lbbe-software.github.io/fitdistrplus/articles/FAQ.html" target="_blank">FAQ</a> before contacting authors

- Marie-Laure Delignette-Muller: `marielaure.delignettemuller<<@))vetagro-sup.fr`
- Christophe Dutang: `dutangc<<@))gmail.com`
- Aurélie Siberchicot: `aurelie.siberchicot<<@))univ-lyon1.fr`

Issues can be reported on [fitdistrplus-issues](https://github.com/lbbe-software/fitdistrplus/issues).

## Citation

If you use `fitdistrplus`, you should cite: <br />
Marie Laure Delignette-Muller, Christophe Dutang (2015). 
*fitdistrplus: An R Package for Fitting Distributions.*
Journal of Statistical Software.
<a href="https://www.jstatsoft.org/article/view/v064i04" target="_blank">https://www.jstatsoft.org/article/view/v064i04</a>
DOI [10.18637/jss.v064.i04](https://doi.org/10.18637/jss.v064.i04).
