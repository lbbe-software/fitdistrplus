
require(actuar)
require(fitdistrplus)
n <- 1e4

truescale <- 100
trueshape1 <- pi

if(FALSE)
{
  plot(ecdf(rburr(n, trueshape1, 2, scale=truescale)))
  lines(ecdf(1/rinvburr(n, trueshape1, 2, scale=1/truescale)), col="red")
  
  curve(pburr(x, trueshape1, 2, scale=truescale), from=0, to=5*truescale)
  curve(pinvburr(1/x, trueshape1, 2, scale=1/truescale, lower=FALSE), add=TRUE, col="red")
}

z <- rburr(n, trueshape1, 2, scale=truescale)
y <- 1/z
fit_y_IB <- fitdist(y, "invburr")
mytitle <- paste("fitted Inv. Burr", paste(signif(coef(fit_y_IB), 3), collapse= ", "))
myleg <- paste("theo. Inv. Burr", paste(signif(c(trueshape1, 2, 1/truescale), 3), collapse= ", "))
cdfcomp(fit_y_IB, xlogscale = TRUE, do.points = FALSE, addlegend = FALSE, main=mytitle)
curve(pinvburr(x, trueshape1, 2, scale=1/truescale), add=TRUE, col="green")
legend("topleft", lty=1, col=c("black", "red", "green"), 
       c("empirical", "fitted", myleg))


y <- rinvburr(n, trueshape1, 2, scale=1/truescale)
fit_y_IB <- fitdist(y, "invburr")
mytitle <- paste("fitted Inv. Burr", paste(signif(coef(fit_y_IB), 3), collapse= ", "))
myleg <- paste("theo. Inv. Burr", paste(signif(c(trueshape1, 2, 1/truescale), 3), collapse= ", "))
cdfcomp(fit_y_IB, xlogscale = TRUE, do.points = FALSE, addlegend = FALSE, main=mytitle)
curve(pinvburr(x, trueshape1, 2, scale=1/truescale), add=TRUE, col="green")
legend("topleft", lty=1, col=c("black", "red", "green"), 
       c("empirical", "fitted", myleg))




