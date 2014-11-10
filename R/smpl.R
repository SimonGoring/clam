#' sample point age estimates from the calibrated distributions ('its' times)
#' 
#' sample point age estimates from the calibrated distributions ('its' times)
#' the probability of a year being sampled is proportional to its calibrated probability
#' @param its
#' @param depths
#' @param calibs
#' @param Est
#' 
#' @return Something
#' @export
#' 
.smpl <- function(its, depths, calibs, Est)
{
  smp <- array(1, dim=c(length(depths), 1+its, 2))
  smp[,1,1] <- Est
  for(i in 1:length(calibs))
    smp[i,(1:its)+1,] <-
    calibs[[i]][sample(1:length(calibs[[i]][,1]), its, prob=calibs[[i]][,2], TRUE),]
  smp
}
