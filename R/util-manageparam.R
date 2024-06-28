#############################################################################
#   Copyright (c) 2024 Christophe Dutang, Marie Laure Delignette-Muller                                                                                                  
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
# checkparam function checks start.arg and fix.arg that parameters are named correctly

# INPUTS 
# start.arg : starting values for optimization or the function to compute them from data or NULL
# fix.arg : fixed values of paramaters or the function to compute them from data or NULL
# obs : the full dataset
# distname : name of the distribution


# OUTPUTS 
# two named list with untested components
manageparam <- function(start.arg, fix.arg, obs, distname)
{
  ddistname <- paste("d", distname, sep="")
  argddistname <- names(formals(ddistname))
  
  #if clause with 3 different cases:
  #start.arg : NULL | named list | a function
  
  if(is.null(start.arg))
  {
    trystart <- try(startargdefault(obs, distname), silent = TRUE)
    if(inherits(trystart, "try-error"))
    {
      cat("Error in computing default starting values.\n")
      stop(trystart)
    }
    lstart <- trystart
    #lstart should be a named list but check it
    hasnoname <- is.null(names(lstart)) || !is.list(lstart)
    if(hasnoname)
      stop("Starting values must be a named list, error in default starting value.")
    
  }else if(is.list(start.arg))
  {
    hasnoname <- is.null(names(start.arg))
    if(hasnoname)
      stop("Starting values must be a named list (or a function returning a named list).")
    lstart <- start.arg
  }else if(is.function(start.arg))
  {
    trystart <- try(start.arg(obs), silent = TRUE)
    if(inherits(trystart, "try-error"))
    {
      cat("Error in computing starting values with your function.\n")
      stop(trystart)
    }
    lstart <- trystart
    hasnoname <- is.null(names(lstart)) || !is.list(lstart)
    if(hasnoname)
      stop("Starting values must be a named list, your function does not return that.")
  }else
    stop("Wrong type of argument for start")
  
  #if clause with 3 different cases:
  #fix.arg : NULL | named list | a function
  if(is.null(fix.arg))
  {
    lfix <- NULL
  }else if(is.list(fix.arg))
  {
    hasnoname <- is.null(names(fix.arg))
    if(hasnoname)
      stop("Fixed parameter values must be a named list (or a function returning a named list).")
    lfix <- fix.arg
  }else if(is.function(fix.arg))
  {
    tryfix <- try(fix.arg(obs), silent = TRUE)
    if(inherits(tryfix, "try-error"))
    {
      cat("Error in computing fixed parameter values with your function.\n")
      stop(tryfix)
    }
    lfix <- tryfix
    hasnoname <- is.null(names(lfix)) || !is.list(lfix)
    if(hasnoname)
      stop("Fixed parameter values must be a named list, your function does not return that.")
  }else
    stop("Wrong type of argument for fix.arg")
  
  #eliminate arguments both in lstart and lfix (when start.arg was NULL)
  if(is.null(start.arg) && !is.null(lfix))
  {
    lstart <- lstart[!names(lstart) %in% names(lfix)]
    if(length(lstart) == 0)
      stop("Don't need to use fitdist() if all parameters have fixed values")
  }
  
  #check if distname has both rate and scale parameter
  if("rate" %in% argddistname && "scale" %in% argddistname)
  {
    if("rate" %in% names(lfix) && "scale" %in% names(lstart))
    {
      lstart <- lstart[names(lstart) != "scale"]
    }else if("scale" %in% names(lfix) && "rate" %in% names(lstart))
    {
      lstart <- lstart[names(lstart) != "rate"]
    }
    if(length(lstart) == 0)
      stop("Don't need to use fitdist() if all parameters have fixed values")
  }
  
  list("start.arg"=lstart, "fix.arg"=lfix)
}

