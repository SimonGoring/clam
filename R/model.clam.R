#' Calculate the age-depth model and its uncertainty
#' 
#' The workhorse function for clam, internal to the package.
#' 
#' @param type numeric.  Model type to generate the chronology.  Default is 1 (linear), also available are 2, polynomial; 3, cubic spline; 4, smooth spline; 5, locally weighted spline (loess).
#' @param smooth
#' @param its
#' @param wghts
#' @param depths
#' @param errors
#' @param depthseq
#' @param prob
#' @param est
#' @param dat
#' @param smp
#' @param greyscale
#' @param remove.reverse
#' @param storedat
#' @param ageofdepth
#' @param BCAD
#' 
#' @return A table of age, depth and error for the core based on a pre-defined method.
#' @export
#' 
.model.clam <- function(type, smooth, its, wghts, depths, errors, depthseq, prob, est, dat, smp, greyscale, remove.reverse, storedat, ageofdepth, BCAD)
{
  # warn for extrapolation, refuse to do so for loess
  if(min(depthseq) < min(dat$depth) || max(depthseq) > max(dat$depth))
    if(type==5)
      stop(" cannot extrapolate using loess! Change settings.\n ", call.=FALSE) else
        cat(" extrapolating beyond dated levels, dangerous!\n ")
  
  # choose model: interpolation, (polynomial) regression, spline, smooth spline or loess
  chron <- array(0, dim=c(length(depthseq), its))
  
  if(type==1)   chron <- .interp(depthseq, depths, its, chron, smp) else
    if(type==2) chron <- .poly(depthseq, smooth, wghts, errors, depths, its, chron, smp) else
    if(type==3) chron <- .spline(depthseq, smooth, depths, its, chron, smp) else
    if(type==4) chron <- .smooth(depthseq, smooth, wghts, errors, depths, its, chron, smp) else
    if(type==5) chron <- .loess(depthseq, smooth, wghts, errors, depths,  its, chron, smp) else
    if(type==6) chron <- .stinterpx(depthseq, smooth, wghts, errors, depths,  its, chron, smp)
  
  # test against age reversals
  warp <- c()
  if(remove.reverse!=FALSE)
    for(i in 1:ncol(chron))
      if(!BCAD && min(diff(chron[,i]), na.rm=TRUE) <= 0 || BCAD && max(diff(chron[,i]), na.rm=TRUE) >= 0)
        warp <- c(warp, i)
  if(length(warp) > 0)
    if(length(warp) > remove.reverse*its)
      cat("\n\n !!! Too many models with age reversals!!!\n") else
      {
        cat("\n Removing", length(warp), "models with age reversals,", its-length(warp), "models left...")
        chron <- chron[,-warp]
        smp <- smp[,-warp,]
      }
  
  if(length(ageofdepth) > 0)
    if(ageofdepth %in% depthseq)
      .ageofdepth <<- chron[which(depthseq==ageofdepth),]
  
  if(storedat)
  {
    chron <<- chron
    smp <<- smp
  }
  
  # find uncertainty ranges of calendar age for each depth of the core
  calrange <- array(0, dim=c(nrow(chron), 2))
  wm <- c()
  for(i in 1:nrow(chron))
  {
    x <- chron[i,2:ncol(chron)]
    qp <- (1-prob)/2
    calrange[i,] <- quantile(x, c(qp, 1-qp))
    if(est==1) wm[i] <- weighted.mean(x)
  }
  if(est==1)
    cbind(depthseq, cbind(calrange, wm)) else
      cbind(depthseq, cbind(calrange, chron[,1]))
}
