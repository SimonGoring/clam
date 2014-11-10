#' Plot age model with proxy data.
#' 
#' Only works after doing a clam run with \code{proxies=TRUE}.
#' 
#' @param prox
#' @param errors=TRUE
#' @param proxcol=grey(0.5)
#' @param revyr=TRUE
#' 
#' @return A plot of the age model function with proxies.
#' @export
#' 
plot.proxies <- function(prox, errors=TRUE, proxcol=grey(0.5), revyr=TRUE)
{
  prx <- dat$proxies
  if(length(prox)>1) layout(matrix(1:length(prox), ncol=1))
  for(j in 1:length(prox))
  {
    pr <- prx[which(!is.na(prx[,prox+1])),]
    ages <- array(0, dim=c(nrow(pr),3))
    for(i in 1:nrow(pr))
      ages[i,] <- calrange[which(calrange[,1]==pr[i,1]),c(2,3,4)]
    xlim <- range(ages)
    if(!dat$BCAD) xlim <- rev(xlim)
    if(revyr) xlim <- rev(xlim)
    plot(ages[,3], pr[,prox+1], type="n", xlim=xlim, xlab=ifelse(dat$BCAD, "cal BC/AD", "cal BP"), ylab=names(pr)[prox+1])
    if(errors)
      for(i in 2:nrow(pr))
        polygon(c(ages[(i-1):i,1], ages[i:(i-1),2]), c(pr[c((i-1):i, i:(i-1)),prox+1]), col=proxcol, border=proxcol)
    lines(ages[,3], pr[,prox+1])
  }
  layout(1)
}
