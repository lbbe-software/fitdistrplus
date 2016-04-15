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
### Infinitely differentiable transformations of R into a bounded or half-bounded interval
###
###         R functions
### 


#inverse are useless?


#Transformation from (-Inf, +Inf) to (-1, 0)
Tm10 <- function(x)
  -1/(1+exp(-x))
#Inverse
iTm10 <- function(x) 
  log(-x/(1+x))

#Transformation from (-Inf, +Inf) to (0, 1)
T01 <- function(x)
  1/(1+exp(-x))
#Inverse
iT01 <- function(x)
  log(x/(1-x))

#Transformation from (-Inf, +Inf) to (0, +Inf)
T0Inf <- function(x)
  exp(x)
#Inverse
iT0Inf <- function(x)
  log(x)

#Transformation from (-Inf, +Inf) to (1, +Inf)
T1Inf <- function(x)
  1+exp(x)
#Inverse
iT1Inf <- function(x)
  log(x-1)

