visualize <- FALSE # TRUE for manual tests with visualization of results

if(visualize)
{
  
  library(parallel)
  
  
  #fonction basique evaluant la moyenne empirique d'un echantillon gaussien
  f <- function(i) mean(rnorm(1e4, mean=i))
  
  #exemple simple
  glist <- 1:20 	
  cl <- parallel::makeCluster(2)
  system.time(
    res <- parallel::parLapply(cl, glist, f)
  )
  parallel::stopCluster(cl)
  
  
  #exemple en faisant varier le nombre de coeurs et le nombre de simulations
  nbsimu <- 10^(1:2)
  cores <- 1:4
  cores <- 1:getOption("cl.cores", 2)
  
  partime <- matrix(NA, length(nbsimu), length(cores)+1)
  colnames(partime) <- c("R", paste("core",cores))
  rownames(partime) <- paste("n", nbsimu, sep="=")
  
  partime[, 1] <- sapply(1:length(nbsimu), function(i) 
    system.time(lapply(1:nbsimu[i], f))[3])
  for(j in 1:length(cores))
  {
    print(cores[j])
    cl <- parallel::makeCluster(cores[j])
    partime[, j+1] <- sapply(1:length(nbsimu), function(i) 
      system.time(parallel::parLapply(cl, 1:nbsimu[i], f))[3])
    parallel::stopCluster(cl)
  }		
  
  
  partime
  
}