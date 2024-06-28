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

#check for uncensored data the presence of NaN, +/-Inf, NA
checkUncensoredNAInfNan <- function(x)
{
  if(!is.numeric(x))
    stop("data is not a numeric vector.")
  if(any(is.nan(x)))
    stop("data contain NaN (not a numeric) values.")
  if(any(is.infinite(x)))
    stop("data contain Inf (infinite) values.")
  if(any(is.na(x)))
    stop("data contain NA values.")
  invisible(NULL)
}

#check for censored data the presence of NaN, +/-Inf, NA
checkCensoredDataFrameNAInfNan <- function(x)
{
  if(!is.data.frame(x))
    stop("censdata is not a dataframe with two columns named left and right.")
  if(NCOL(x) != 2)
    stop("censdata is not a dataframe with two columns named left and right.")
  if(!"left" %in% colnames(x) || !"right" %in% colnames(x))
    stop("censdata is not a dataframe with two columns named left and right.")
  if(any(!is.numeric(x$left) | !is.numeric(x$right)))
    stop("censdata contain NaN (not a numeric) values.")
  if(any(is.nan(x$left) | is.nan(x$right)))
    stop("censdata contain NaN (not a numeric) values.")
  if(any(is.infinite(x$left) | is.infinite(x$right)))
    stop("censdata contain Inf (infinite) values.")
  if(any(is.na(x$left) & is.na(x$right)))
    stop("censdata contain two NA values on the same line.")
  leftsupright <- x$left > x$right
  leftsupright <- leftsupright[!is.na(leftsupright)]
  if (any(leftsupright))
    stop("each censdata$left value must be less or equal to the corresponding censdata$right value")
  invisible(NULL)
}
