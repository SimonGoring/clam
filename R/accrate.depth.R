#' Accumulation rate calculation and plotting.
#' 
#' Calculate the slope of a straight curve between depths above and below the desired point (for each iteration). Requires sufficiently dense density of depths, e.g. yrsteps=1
#' to calculate accumulation rates at a depth. Before running this, run your core in clam and store the data, so, provide the option storedat=TRUE#' This function is internal to Clam.
#' 
#' @param depth numeric.  The radiocarbon year.
#' @param yrcm boolean. The estimated error for the <sup>14</sup>C age.
#' @param prob numeric. The confidence level (0 - 1).
#' 
#' @return A plot of the accumulation rate by depth.
#' @export

accrate.depth <- function(depth, yrcm=TRUE, prob=.95)
{
  if(depth <= min(calrange[,1]) || depth >= max(calrange))
    stop("Accumulation rates cannot be calculated for the top or bottom of the core. Please check the manual", call.=FALSE)
  d <- max(which(calrange[,1] <= depth))
  if(yrcm)
    accrate <- (chron[d+1,]-chron[d-1,]) / (calrange[d+1,1]-calrange[d-1,1]) else
      accrate <- (calrange[d+1,1]-calrange[d-1,1]) / (chron[d+1,]-chron[d-1,])
  acc <- density(accrate)
  plot(acc, main="", xlab=if(yrcm) "yr/cm" else "cm/yr")
  abline(h=0)
  o <- order(acc$y, decreasing=TRUE)
  acc <- cbind(acc$x[o], cumsum(acc$y[o])/sum(acc$y))
  acc <- range(acc[acc[,2] <= prob,1])
  rect(acc[1], 0, acc[2], -999, col=grey(.5), border=grey(.5))
  cat(100*prob, "% ranges: ", acc[1], " to ", acc[2], if(yrcm) " yr/cm\n" else " cm/yr\n", sep="")
}
