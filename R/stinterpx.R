#' Smoothing using Stineman interpolation.
#' 
#' Fitting a cubic smoothed splines through the data, with smoothing factor.
#' 
#' @param depthseq 
#' @param smooth
#' @param wghts
#' @param errors
#' @param depths
#' @param its
#' @param chron
#' @param smp
#' 
#' @return The smoothed chronology.
#' @export
#' 

.stinterpx <- function(depthseq, smooth, wghts, errors, depths, its, chron, smp)
{
    
    require(stinepack)
    
    for(i in 1:its){
      if(wghts==1) w <- smp[,i,2]
      
      chron[,i] <- stinterp(x = depths, y = smp[,i,1], xout = depthseq)  
      if(i/(its/5) == round(i/(its/5))) cat(".")
    }
    
    chron
    
  }
