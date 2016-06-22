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
### extract model coefficients 
###
###         R functions
###

#already in R
#
#coef <- function(object, ...)
#    UseMethod("coef")
#
#coef.default <- function(object, ...)
#    return(object)

coef.fitdist <- function(object, ...)
{
    stopifnot(inherits(object, "fitdist"))  
    
    if(is.null(object$estimate))
        stop("Internal error in coef.fitdist")
    else
        return(object$estimate)
}

coef.fitdistcens <- function(object, ...)
{
    stopifnot(inherits(object, "fitdistcens"))
    
    if(is.null(object$estimate))
        stop("Internal error in coef.fitdistcens")
    else
        return(object$estimate)
}
