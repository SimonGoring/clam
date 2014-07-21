#' Interpolate linearly between the data (default)
#' 
#' Provides linear interpolation at all depths between the \code{depthseq} over a large set of runs, defined by \code{its}.
#' 
#' @param depthseq
#' @param depths
#' @param its
#' @param chron
#' @param smp
#' 
#' @return A vector of dates.
#' @export


.interp <- function(depthseq, depths, its, chron, smp)
{
  cat(" Interpolating, sampling")
  for(i in 1:its)
  {
    temp <- approx(depths, smp[,i,1], depthseq, ties=mean)$y
    
    # allow for extrapolation... dangerous!
    if(min(depthseq) < min(depths))
    {
      minus <- which(depthseq < min(depths))
      slope <- diff(temp)[max(minus)+1]/diff(depthseq)[max(minus)+1]
      temp[minus] <- temp[max(minus)+1] + slope * (depthseq[minus] - min(depths))
    }
    if(max(depthseq) > max(depths))
    {
      maxim <- which(depthseq > max(depths))
      slope <- diff(temp)[min(maxim)-2]/diff(depthseq)[min(maxim)-2]
      temp[maxim] <- temp[min(maxim)-1] + slope * (depthseq[maxim] - max(depths))
    }
    chron[,i] <- temp
    if(i/(its/5) == round(i/(its/5))) cat(".")
  }
  chron
}
