#' Calculate pMC values from 14C ages.
#' 
#' Calculate the slope of a straight curve between depths above and below the desired point (for each iteration). Requires sufficiently dense density of depths, e.g. yrsteps=1
#' to calculate accumulation rates at a depth. Before running this, run your core in clam and store the data, so, provide the option storedat=TRUE#' This function is internal to Clam.
#' 
#' @param mn numeric.
#' @param sdev numeric.
#' @param ratio numeric.
#' @param decimals numeric.
#' 
#' @return Rounded values for ???mn??? to a specified decimal level.
#' @export

age.pMC <- function(mn, sdev, ratio=100, decimals=3)
{
  y <- exp(-mn/8033)
  sdev <- y - exp(-(mn+sdev)/8033)
  signif(ratio*c(y, sdev), decimals)
}
