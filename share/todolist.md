*TODO list*

# Shorter term 


1. [ ] take into account weights for definition of starting values
1. [ ] take into account weights in every functions and add examples in fitdist… gofstat, plotdist, plotdistcens, descdist, and all plotting functions
1. [ ] ajouter des choix de valeurs initiales pour des lois de actuar, vgam (packages considérés comme centro dans la task view (mention core dans lis à la fin) et bien utilisés avec le format classique, actuar, vgam, gamlss.dist) et faire une FAQ associée avec un tableau lisant les dist prises en compte et la méthode associées (moments ou quantiles), voire la formule
1. [ ] Traiter dans la FAQ la question du choix du nombre d'itérations bootstrap - compléter avec un 4.4 en donnant un exemple où on fait varier le nb d'itérations (faut que ça se stabilise)
1. [ ] Add stats for fits on censored data and the corresponding  gofstat function : look at recent papers
1. [ ] MSE for censored data ?
1. [ ] explore Cullen and Frey for various dist with trimmed linear moments
1. [ ] look in scholar google : fitdistrplus "survival data" and fitdistrplus dependencies to appreciate the use of fitdistrplus on survival data and the needs utiliser le package Rwsearch
1. [ ] biblio sur critères d'information (AICc dans le cas général ou autres critères) pour voir si on peut élargir les stat données
1. [ ] faire un script automatisé pour analyser les résumés et les revues des articles citant fitdistrplus
1. [ ] add calculation of the hession using optimHess within fitdist or mledist in cases where it is not given by optim
1. [ ] better compute the hessian matrix for MLE (add one step after the final estimate to compute correctly hessian, see example for gamma distribution)
1. [ ] think about a hexa-logo (histogram with fitted densities on a simulated example from two asymetric  distribution (lnorm, weibull) without axis)
1. [ ] add a cheatsheet (refcard) for fitdistrplus (possibly in pptx?), think about sections
1. [ ] make a flow chart on how functions are interrelated (see techdoc/)
1. [X] Make a markdown TODO list
1. [ ] delete the R forge project https://r-forge.r-project.org/projects/riskassessment/ ?
1. [ ] should we return the points/lines/rectanges drawn in invisible()? as in hist.default() called by hist() ou plot.stepfun called by plot(ecdf())? 


# Of less priority

1. [ ] add the Wasserstein-Kantorovich distance in mge and gofstat (see Del Barrio 1999 and 2000 and Gibbs 2002)
1. [ ] offer a function to do prior elicitation from quantiles (look before at package expert)
1. [ ] multivariate distribution fitting, in particular copula fitting
1. [ ] build a function to fit a same distribution to a variable (coded on one column) on different subsets defined by a factor (coded on another column) : tools for exploration before modeling the effect of the factor on distribution parameters (dans un autre package ? est-ce déjà fait ailleurs?) gfitdistrplus : a new package to build
1. [ ] bibliométrie de l'utilisation du package dans les publis récentes (données censurées, graphes d'ajustement….) depuis la soumission JSS


# Longer term 
