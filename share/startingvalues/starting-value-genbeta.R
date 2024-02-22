library(actuar)

n <- 1e3
x <- rgenbeta(n, 2, 2, 2, scale=2)
x <- rgenbeta(n, 2, 2, 2, scale=1/2)

f1 <- function(x, a, b)
{
  gamma(a+1/x)/gamma(a+b+1/x)
}

curve(f1(x, 1, 1), from=.5, to =10)

theta <- mean(x^2)/mean(x)
y <- (x*theta/max(x))^3

cdfcomp(fitdist(y, "beta"))
