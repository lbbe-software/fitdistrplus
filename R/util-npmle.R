# ----------------------------------------------------------------------- #
#   Copyright (c) 2020 Marie Laure Delignette-Muller                      #
#                      and Christophe Dutang
#                                                                         #
# Nonparametric maximum likelihood estimation from interval-censored data #
# ----------------------------------------------------------------------- #
# Function calling alternatives to npsurv() from the npsurv package       #
# ----------------------------------------------------------------------- #

  
npmle <- function(censdata, method = "Wang")
{
  method <- match.arg(method, c("Wang", "Turnbull.intervals"))
  if (method == "Wang")
  {
    db <- censdata
    db$left[is.na(db$left)] <- -Inf
    db$right[is.na(db$right)] <- Inf
    
    ############## TO CHANGE IF npsurv.minimal call is changed
    r <- npsurv.minimal(db, pkg="stats")
    if (r$convergence)
    {
      f <- r$f
    } else
    {
      warning("Due to lack of convergence of Wang algorithm, method Turnbull.intervals 
              was used instead for NPMLE")
      f <- Turnbull.intervals(censdata)
    }
  } else
  if (method == "Turnbull.intervals")
  {
    f <- Turnbull.intervals(censdata)
  }
  return(f)
  
}
