#' @title Accumulation rate calculation and plotting.
#' 
#' Calculate the accumulation rate (slope of a straight curve between depths) above and below the desired point. Requires sufficiently dense density of depths, e.g. steps=1
#' to calculate accumulation rate at an age. Before doing this, run your core in clam and store the data, so, provide the option storedat=TRUE
#' This function is internal to Clam.
#' 
#' @param age numeric.  The calibrated year.
#' @param yrcm boolean. The estimated error for the <sup>14</sup>C age.
#' @param prob numeric. The confidence level (0 - 1).
#' 
#' @return A plot of the accumulation rate by age.
#' @export

accrate.age <- function(age, yrcm=TRUE, prob=.95)
{
  
  accrate <- rep(NA, ncol(chron))
  
  for(i in 1:ncol(chron))
  {
    a <- max(which(chron[,i] <= age))
    if(yrcm){
      accrate[i] <- (chron[a+1,i]-chron[a-1,i]) / (calrange[a+1,1]-calrange[a-1,1]) 
    }
    else {
      accrate[i] <- (calrange[a+1,1]-calrange[a-1,1]) / (chron[a+1,i]-chron[a-1,i])
    }
  }

  acc <- density(accrate)
  
  plot(acc, main="", xlab=if(yrcm) "yr/cm" else "cm/yr")
  abline(h=0)
  
  o <- order(acc$y, decreasing=TRUE)
  acc <- cbind(acc$x[o], cumsum(acc$y[o])/sum(acc$y))
  acc <- range(acc[acc[,2] <= prob,1])
  
  rect(acc[1], 0, acc[2], -999, col=grey(.5), border=grey(.5))
  cat(100*prob, "% ranges: ", acc[1], " to ", acc[2], if(yrcm) " yr/cm\n" else " cm/yr\n", sep="")  

}