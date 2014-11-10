#' Smoothing cubic splines for the age model.
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
.smooth <- function(depthseq, smooth, wghts, errors, depths, its, chron, smp)
{
  if(length(smooth) < 1) smooth <- .3
  cat(paste(" Using smoothing spline (smoothing ", smooth, "), sampling", sep=""))
  if(wghts==0) w <- c() else w <- 1/errors^2
  for(i in 1:its)
  {
    if(wghts==1) w <- smp[,i,2]
    chron[,i] <- predict(smooth.spline(depths, smp[,i,1], w=w, spar=smooth), depthseq)$y
    if(i/(its/5) == round(i/(its/5))) cat(".")
  }
  chron
}