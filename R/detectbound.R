#############################################################################
#   Copyright (c) 2016 Marie Laure Delignette-Muller, Christophe Dutang                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################
### detect boundaries

detectbound <- function(distname, vstart, obs, fix.arg=NULL, echo=FALSE)
{
  ddistname <- paste("d", distname, sep="")
  argdist <- formalArgs(ddistname)
  argdist <- argdist[!argdist %in% c("x", "log")]
  
  stopifnot(all(names(vstart) %in% argdist))
  
  if("scale" %in% argdist && "rate" %in% argdist)
  {
    if(length(grep("rate", as.character(formals(ddistname)$scale))) > 0)
    {
      argdist <- argdist[argdist != "rate"] #remove rate for parameter list
      if("rate" %in% names(vstart)) #update value if rate is supplied
      {
        vstart["rate"] <- 1/vstart["rate"]
        names(vstart)[names(vstart) == "rate"] <- "scale"
      }
    }
  }
  argdist <- argdist[!argdist %in% names(fix.arg)]
  
  if(length(argdist) == 0)
    return(NULL)
  
  if(echo)
  {
    print(argdist)
    print(vstart)
  }
  lowb <- rep(-Inf, length(argdist))
  uppb <- -lowb
  names(lowb) <- names(uppb) <- argdist
  
  eps <- sqrt(.Machine$double.eps)
  
  owarn <- getOption("warn")
  oerr <- getOption("show.error.messages")
  options(warn=-1, show.error.messages=FALSE)
  
  for(a in argdist)
  {
    if(echo)
      cat(a, "\n")
    dx <- do.call(ddistname, c(list(obs), as.list(vstart), as.list(fix.arg)))
    if(any(is.nan(dx)))
      stop("wrong init param")
    vstarta <- vstart
    aval <- -1:1
    
    for(i in 1:length(aval))
    {
      vstarta[a] <- aval[i]-eps
      dx1 <- try(do.call(ddistname, c(list(obs), as.list(vstarta), as.list(fix.arg))), silent=TRUE)
      vstarta[a] <- aval[i]+eps
      dx2 <- try(do.call(ddistname, c(list(obs), as.list(vstarta), as.list(fix.arg))), silent=TRUE)
      if(echo)
      {
        print(dx1)
        print(dx2)
      }
      
      if(class(dx1) == "try-error" && class(dx2) != "try-error")
      {
        lowb[a] <- aval[i]
      }
      if(all(is.nan(dx1)) && all(!is.nan(dx2)))
      {
        lowb[a] <- aval[i]
      }
      if(class(dx1) != "try-error" && class(dx2) == "try-error")
      {
        uppb[a] <- aval[i]
      }
      if(all(!is.nan(dx1)) && all(is.nan(dx2)))
      {
        uppb[a] <- aval[i]
      }
    }
      
  }
  
  options(warn=owarn, show.error.messages=oerr)
  rbind(lowb, uppb)
}

