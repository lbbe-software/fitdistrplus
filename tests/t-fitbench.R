require(fitdistrplus)


if(FALSE)
{  
  
  
  x <- rgamma(1e3, shape=3/2, rate= 1/2)
  
  initval <- unlist(fitdistrplus:::start.arg.default(x, "gamma"))
  
  fitdistrplus:::fitbench(x, "gamma", "mle")
  
  fitdistrplus:::fitbench(x, "gamma", "mle", lower=0)
  
  
  # grgam(1:2, x)
  # grgamlnl(x[1], shape=1, rate=2)
  # grgamlnl(x[1:2], shape=1, rate=2)
  
  fitdistrplus:::fitbench(x, "gamma", "mle", grad=fitdistrplus:::grlnlgamma)
  
  fitdistrplus:::fitbench(x, "gamma", "mle", grad=fitdistrplus:::grlnlgamma, lower=0)
  
  #mledist(x, "gamma", grad=grgam, lower=0, optim.method = "CG")
  #mledist(x, "gamma", grad=grgam, lower=0, optim.method = "BFGS")
  
}
