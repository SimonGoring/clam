#' LOESS Interpolation at depth.
#' 
#' Fits a locally weighted (1/errors^2) splines through the data, with smoothing factor
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
#' @return A vector of dates.
#' @export
.loess <- function(depthseq, smooth, wghts, errors, depths, its, chron, smp)
{
  if(length(smooth) < 1) smooth <- .75
  cat(paste(" Using loess (smoothing ", smooth, "), sampling", sep=""))
  if(wghts==0) w <- c() else w <- 1/errors^2
  for(i in 1:its)
  {
    if(wghts==1) w <- smp[,i,2]
    chron[,i] <- predict(loess(smp[,i,1] ~ depths, weights=w, span=smooth), depthseq)
    if(i/(its/5) == round(i/(its/5))) cat(".")
  }
  chron
}