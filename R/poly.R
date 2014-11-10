#' Polynomial age-depth regression.
#' 
#' A polynomial regressions of certain order through the data (default linear, y=ax+b)
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
#' @return A chronology based on polynomial regression.
#' @export
#' 
.poly <- function(depthseq, smooth, wghts, errors, depths, its, chron, smp)
{
  if(length(smooth)==0)
    cat(" Using linear regression, sampling") else
      cat(paste(" Using polynomial regression (degree ", smooth, "), sampling", sep=""))
  if(wghts==0) w <- c() else w <- 1/errors^2
  for(i in 1:its)
  {
    if(wghts==1) w <- smp[,i,2]
    chron[,i] <- predict(lm(smp[,i,1] ~ poly(depths, max(1, smooth)), weights=w), data.frame(depths=depthseq))
    if(i/(its/5) == round(i/(its/5))) cat(".")
  }
  chron
}
