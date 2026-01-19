# *TODO list*

## Shorter term

1.  [ ] take into account weights for definition of starting values
2.  [ ] take into account weights in every functions and add examples in fitdist, gofstat, plotdist, plotdistcens, descdist, and all plotting functions
3.  [x] ajouter des choix de valeurs initiales pour des lois de actuar, vgam (packages considérés comme centro dans la task view (mention core dans lis à la fin) et bien utilisés avec le format classique, actuar, vgam, gamlss.dist) et faire une FAQ associée avec un tableau lisant les dist prises en compte et la méthode associées (moments ou quantiles), voire la formule
4.  [x] Traiter dans la FAQ la question du choix du nombre d'itérations bootstrap - compléter avec un 4.4 en donnant un exemple où on fait varier le nb d'itérations (faut que ça se stabilise)
5.  [ ] Add stats for fits on censored data and the corresponding gofstat function : look at recent papers
6.  [ ] consider the method MSE for censored data ?
7.  [ ] explore Cullen and Frey for various dist with trimmed linear moments as well as other common distributions
8.  [ ] look in scholar google : fitdistrplus "survival data" and fitdistrplus dependencies to appreciate the use of fitdistrplus on survival data and the needs utiliser le package Rwsearch
9.  [ ] biblio sur critères d'information (AICc dans le cas général ou autres critères) pour voir si on peut élargir les stat données
10. [ ] faire un script automatisé pour analyser les résumés et les revues des articles citant fitdistrplus
11. [x] add calculation of the hessian using `optimHess` within `fitdist` or `mledist` in cases where it is not given by optim *CD*
12. [x] better compute the hessian matrix for MLE (add one step after the final estimate to compute correctly hessian, see example for gamma distribution)
13. [x] think about a hexa-logo (histogram with fitted densities on a simulated example from two asymetric distribution (lnorm, weibull) without axis)
14. [ ] add a cheatsheet (refcard) for fitdistrplus (possibly in pptx?), think about sections, update the website github.io for the Reference page
15. [ ] make a flow chart on how functions are interrelated (see techdoc/)
16. [x] Make a markdown TODO list
17. [ ] delete the R forge project <https://r-forge.r-project.org/projects/riskassessment/> ?
18. [x] should we return the points/lines/rectangles drawn in `invisible()`? as in `hist.default()` called by `hist()` ou `plot.stepfun` called by `plot(ecdf())`? : only when `plotstyle == "graphics"`
19. [x] add a new question to the FAQ, following <https://github.com/lbbe-software/fitdistrplus/issues/21>
20. [x] add a new question to the FAQ, following <https://github.com/lbbe-software/fitdistrplus/issues/23>
21. [x] add a generic function `AIC` and `BIC` for `fitdist(cens)` objects
22. [x] close following issues #2, #21, #23, #25 appropriately
23. [x] improve Cullen-Frey graph following <https://github.com/lbbe-software/fitdistrplus/issues/24>
24. [x] manage differential graphical parameters in `plotdist` see <https://github.com/lbbe-software/fitdistrplus/issues/27>, link FAQ 5.2
25. [x] make a `gofstat` for `fitdistcens` objects *MLDM* et *CD* : create `gofstat.fitdistcens`
26. [x] make introduction for vignette `starting values` and update layout, update README on `github.io` *CD*
27. [x] update `mledist.Rd` with new starting value list (all distribution in `actuar` except phase-type), link with vignette *CD*
28. [x] check that all base R distribution have starting values *CD*
29. [x] add defensive programming for infinite / NaN / NA values for `fitdist` objects *MLDM* et *CD*
30. [x] add defensive programming for NaN values for `fitdistcens` objects, convert `Inf` to `NA` (raise error for inconsistent values), raise error for double `NA` *MLDM* et *CD* 
31. [x] add a link to `Surv2fitdistcens` in `fitdistcens`, add `Surv2fitdistcens` when infinite values
32. [x] check link of `Surv2fitdistcens` in the FAQ
33. [ ] add exploratory tools for assessing the role of covariates and for choosing the appropriate distribution in modelling (survival parametric model), e.g. dose-response
34. [x] allow the computation of vcov for MME, remove NOTE when checking the package *CD*
35. [x] add DOI for reference in man pages
36. [ ] allow the use of a function for `distr` argument rather than a character string, e.g., `evir::dgev` so that `evir` is not attached. probably designed for functions `mledist`, `mmedist`, `qmedist`,... rather than for `fitdist`. See Issue https://github.com/lbbe-software/fitdistrplus/issues/38
37. [x] move some test files in `share/`, e.g., `t-starting-values*` *CD*
38. [ ] harmonize the use of `set.seed()` in test files and in man files *MLDM*
39. [x] close Issue 38 <https://github.com/lbbe-software/fitdistrplus/issues/38> with a comment.



## Longer term

1.  [ ] add the Wasserstein-Kantorovich distance in `mge` and `gofstat` (see Del Barrio 1999 and 2000 and Gibbs 2002)
2.  [ ] offer a function to do prior elicitation from quantiles (look before at package expert)
3.  [ ] multivariate distribution fitting, in particular copula fitting
4.  [ ] build a function to fit a same distribution to a variable (coded on one column) on different subsets defined by a factor (coded on another column) : tools for exploration before modeling the effect of the factor on distribution parameters (dans un autre package ? est-ce déjà fait ailleurs?) gfitdistrplus : a new package to build
5.  [ ] bibliométrie de l'utilisation du package dans les publis récentes (données censurées, graphes d'ajustement.) depuis la soumission JSS



