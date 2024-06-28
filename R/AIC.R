#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Christophe Dutang                                                                                                  
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
### AIC generic functions
###
###         R functions
###

AIC.fitdist <- function(object, ..., k = 2)
{
  stopifnot(inherits(object, "fitdist"))
  
  if(is.null(object$aic))
    stop("Internal error in AIC.fitdist")
  else
  {
    if (k == 2) 
    {
      return(object$aic)
    } else
    {
      npar <- object$aic / 2 +  object$loglik
      aic <- -2 * object$loglik + k * npar
      return(aic)
    }
  }
}

AIC.fitdistcens <- function(object, ..., k = 2)
{
  stopifnot(inherits(object, "fitdistcens"))
  
  if(is.null(object$aic))
    stop("Internal error in AIC.fitdistcens")
  else
  {
    if (k == 2) 
    {
      return(object$aic)
    } else
    {
      npar <- object$aic / 2 +  object$loglik
      aic <- -2 * object$loglik + k * npar
      return(aic)
    }
  }
}
