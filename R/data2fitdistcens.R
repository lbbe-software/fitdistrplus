#############################################################################
#   Copyright (c) 2021 Christophe Dutang                                                                                                  
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
### data handling functions
###
###         R functions
### 

data2fitdistcens <- function(time, time2, event,
                             type=c('right', 'left', 'interval'))
{
  type <- match.arg(type,  c('right', 'left', 'interval'))
  stopifnot(is.numeric(time))
  stopifnot(is.numeric(time2))
  stopifnot(length(time) >= 1)
  stopifnot(length(time) == length(time2))
  stopifnot(length(time) == length(event))
  
  if(is.logical(event))
    event <- 1*(event == TRUE)
  else if(all(event == 1 | event == 2))
  {
    event <- 1*(event == 2)
    stopifnot(type %in% c('right', 'left'))
  }else if(all(event %in% 0:3) && sum(event >= 3) > 0)
  {
    stopifnot(type %in% 'interval')
  }else if(is.factor(event))
  {
    stopifnot(length(levels(event)) == 2)
    event <- 1*(event == levels(event)[2])
  }
  if(any(!event %in% 0:1) && type != "interval")
    stop("wrong 'event' argument")
  
  #compute data.frame
  if(type == "right")
  {
    out <- cbind.data.frame(left=time,
                            right=NA)
    out$right[event == 1] <- time2[event == 1]
  }else if(type == "left")
  {
    out <- cbind.data.frame(left=NA,
                            right=time2)
    out$left[event == 1] <- time[event == 1]
  }else #type == "interval"
  {
    out <- cbind.data.frame(left=rep(NA, length(time)),
                            right=rep(NA, length(time2)))
    #0=right censored, 
    out$right[event == 0] <- time2[event == 0]
    #1=event at time, 
    out$left[event == 1] <- time[event == 1]
    out$right[event == 1] <- time[event == 1]
    #2=left censored, 
    out$left[event == 2] <- time[event == 2]
    #3=interval censored
    out$left[event == 3] <- time[event == 3]
    out$right[event == 3] <- time2[event == 3]
  }
  out
}
