
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



if(FALSE)
{
z <- rinvburr(n, trueshape1, 2, scale=truescale)
y <- 1/z
fit_y_B <- fitdist(y, "burr", lower=0)
mytitle <- paste("fitted Burr", paste(signif(coef(fit_y_B), 3), collapse= ", "))
myleg <- paste("theo. Burr", paste(signif(c(trueshape1, 2, 1/truescale), 3), collapse= ", "))

cdfcomp(fit_y_B, xlogscale = TRUE, do.points = FALSE, addlegend = FALSE, main=mytitle)
curve(pburr(x, trueshape1, 2, scale=1/truescale), add=TRUE, col="green")
legend("topleft", lty=1, col=c("black", "red", "green"),
       c("empirical", "fitted", myleg))


y <- rburr(n, trueshape1, 2, scale=1/truescale)
fit_y_B <- fitdist(y, "burr", lower=0)
mytitle <- paste("fitted Burr", paste(signif(coef(fit_y_B), 3), collapse= ", "))
myleg <- paste("theo. Burr", paste(signif(c(trueshape1, 2, 1/truescale), 3), collapse= ", "))

cdfcomp(fit_y_B, xlogscale = TRUE, do.points = FALSE, addlegend = FALSE, main=mytitle)
curve(pburr(x, trueshape1, 2, scale=1/truescale), add=TRUE, col="green")
legend("topleft", lty=1, col=c("black", "red", "green"),
       c("empirical", "fitted", myleg))

}






trueshape1 <- 1/2
trueshape2 <- 4
truescale <- 1000

y <- rinvburr(n, trueshape1, trueshape2, scale=truescale)
fitdistrplus:::startarg_fellerpareto_family(y, "invburr")

fit_y_IB <- fitdist(y, "invburr", lower=0)
mystart <- fitdistrplus:::startarg_fellerpareto_family(y, "invburr")
mystart$scale <- 1/mystart$scale
fit_y_IB_newstart1 <- fitdist(y, "invburr", start=mystart, lower=0)
mystart2 <- list("shape1"=1, "shape2"=1, "scale"=1)
fit_y_IB_newstart2 <- fitdist(y, "invburr", start=mystart2, lower=0)
mystart3 <- mystart
mystart3$shape1 <- 10
fit_y_IB_newstart3 <- fitdist(y, "invburr", start=mystart3, lower=0)

mytitle <- paste("fitted Inv. Burr", paste(signif(coef(fit_y_IB), 3), collapse= ", "))
myleg <- paste("theo. Inv. Burr", paste(signif(c(trueshape1, 2, truescale), 3), collapse= ", "))
cdfcomp(list(fit_y_IB, fit_y_IB_newstart2, fit_y_IB_newstart3), xlogscale = TRUE, do.points = FALSE, 
        addlegend = TRUE, main=mytitle, fitlwd=1.5)
curve(pinvburr(x, trueshape1, 2, scale=truescale), add=TRUE, col="blue")
legend("topleft", lty=1, col=c("black", "red", "green", "blue"), 
       c("empirical", "fitted", "fitted", myleg))





x <- rburr(n, trueshape1, trueshape1, scale=1/1e5)
fitdistrplus:::startarg_fellerpareto_family(x, "burr")

fit_y_B <- fitdist(x, "burr", lower=1e-2)

mytitle <- paste("fitted Burr", paste(signif(coef(fit_y_B), 3), collapse= ", "))
myleg <- paste("theo. Burr", paste(signif(c(trueshape1, 2, 1/truescale), 3), collapse= ", "))

cdfcomp(fit_y_B, xlogscale = TRUE, do.points = FALSE, addlegend = FALSE, main=mytitle)
curve(pburr(x, trueshape1, trueshape1, scale=1/1e5), add=TRUE, col="green")
legend("topleft", lty=1, col=c("black", "red", "green"),
       c("empirical", "fitted", myleg))




