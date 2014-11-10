#' Fit cubic spline interpolations through the data
#' 
#' @param depthseq
#' @param smooth
#' @param depths
#' @param its
#' @param chron
#' @param smp
#' 
#' @return A chronology based on cubic spline fitting in the \code{stats} package.
#' @export
#' 
.spline <- function(depthseq, smooth, depths, its, chron, smp)
{
  if(length(smooth) < 1) smooth <- .3
  cat(paste(" Using cubic spline sampling", sep=""))
  for(i in 1:its)
  {
    chron[,i] <- stats::spline(depths, smp[,i,1], xout=depthseq)$y
    if(i/(its/5) == round(i/(its/5))) cat(".")
  }
  chron
}
