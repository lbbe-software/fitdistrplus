#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Regis Pouillot, Jean-Baptiste Denis, Christophe Dutang                                                                                                  
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
### log likelihood
###
###         R functions
###

#logLik already defined in R
#
#logLik <- function(object, ...)
#    UseMethod("logLik")
#
#logLik.default <- function(object, ...)
#    return(object)

logLik.fitdist <- function(object, ...)
{
    stopifnot(class(object) == "fitdist")
    
    if(is.null(object$loglik))
        stop("Internal error in loglik.fitdist")
    else
        return(object$loglik)
}

logLik.fitdistcens <- function(object, ...)
{
    stopifnot(class(object) == "fitdistcens")
    
    if(is.null(object$loglik))
        stop("Internal error in loglik.fitdistcens")
    else
        return(object$loglik)
}
