if(FALSE)
{
  
  require(fitdistrplus)
  
  
  #test actuar initialization
  
  start.arg.default <- fitdistrplus:::start.arg.default
  
  #burr
  library(actuar)
  alpha <- 3
  x <- rburr(1000, alpha, 4, .1)
  
  initburr <- function(x)
  {
    pi <- 1:3/4
    qi <- 1-pi
    xi <- as.numeric(quantile(x, probs=pi))
    y <- log(xi[2])/log(xi[1]/xi[2])
    y1 <- log(xi[3]/xi[2])/log(xi[2]/xi[1])
    y2 <- log(xi[1]/xi[3])/log(xi[2]/xi[1])
    f <- function(eta)
      (qi[1]^eta-1)^y1*(qi[2]^eta-1)^y2*(qi[3]^eta-1) - 1
    eta <- try(uniroot(f, c(-10, -1e-6))$root, silent=TRUE)
    if(class(eta) == "try-error")
      eta <- -1
    alpha <- -1/eta
    lambda <- (qi[1]^eta-1)^y*(qi[2]^eta-1)^(-y-1)
    gamma <- log(lambda*(qi[1]^eta-1))/log(xi[1])
    theta <- lambda^(1/gamma)
    list(shape1=alpha, shape2=gamma, rate=1/theta)
  }
  
  initburr(x)
  
  #fitdist(x, "burr", lower=0, start=initburr(x))
  
  #transformed gamma
  x <- rtrgamma(1000, 3, 2, .1)
  
  start1 <- start.arg.default(x, "trgamma")
  
  plot(ecdf(x))
  curve(ptrgamma(x, start1$shape1, start1$shape2, start1$rate), add=TRUE, lty=2)
  
  #fitdist(x, "trgamma", lower=0)
  
  #inverse transformed gamma
  y <- rinvtrgamma(1000, 3, 2, .1)
  
  start1 <- start.arg.default(y, "invtrgamma")
  
  plot(ecdf(y))
  curve(pinvtrgamma(x, start1$shape1, start1$shape2, start1$rate), add=TRUE, lty=2)
  
  #fitdist(y, "invtrgamma")
  
  
}
