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
#' @importFrom stinepack stinterp
#' @export
#' 

.stinterpx <- function(depthseq, smooth, wghts, errors, depths, its, chron, smp)
{
  
  for(i in 1:its){
      if(wghts==1) w <- smp[,i,2]
      
      chron[,i] <- stinterp(x = depths, y = smp[,i,1], xout = depthseq, method = 'stineman')$y
      
      if(any(depthseq > max(depths))){
        last.good <- data.frame(x = rev(depthseq[!is.na(chron[,i])])[1:2],
                                y = rev(chron[!is.na(chron[,i]),i])[1:2])
        reset <- predict(lm(y ~ x, data = last.good), 
                         newdata = data.frame(x = depthseq[is.na(chron[,i]) & depthseq > max(depths)]))
        chron[is.na(chron[,i]) & depthseq > max(depths),i] <- reset
      }
      
      if(any(depthseq < min(depths))){
        last.good <- data.frame(x = depthseq[!is.na(chron[,i])][1:2],
                                y = chron[!is.na(chron[,i]),i][1:2])
        reset <- predict(lm(y ~ x, data = last.good),
                         newdata = data.frame(x = depthseq[is.na(chron[,i]) & depthseq < min(depths)]))
        chron[is.na(chron[,i]) & depthseq < min(depths),i] <- reset
      }
      if(i/(its/5) == round(i/(its/5))) cat(".")
    }
    
    chron
    
  }
