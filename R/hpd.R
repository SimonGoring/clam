#' Highest posterior density.
#' 
#' Find the highest posterior density (hpd) of the calibrated distribution.
#' @param dat
#' @param prob
#' @param hpdsteps
#' @param yrsteps
#' 
#' @return The value of the highest posterior density.
#' @export

.hpd <- function(dat, prob, hpdsteps, yrsteps)
{
  # interpolate and rank the ages according to their calibrated distribution probabilities
  dat <- approx(dat[,1], dat[,2], seq(min(dat[,1]), max(dat[,1]), by=yrsteps))
  o <- order(dat$y, decreasing=TRUE)
  dat <- cbind(dat$x[o], dat$y[o]/sum(dat$y))
  
  # only retain those ages with cumulative normalised probabilities within required percentage
  dat <- dat[which(cumsum(dat[,2]) <= prob),]
  dat <- dat[order(dat[,1]),]
  
  # identify any individual ranges within the hpd range and calculate their probability
  dif <- which(diff(dat[,1]) > hpdsteps)
  if(length(dif)==0)
    hpds <- cbind(min(dat[,1]), max(dat[,1]), 100*prob) else
    {
      dif <- c(dat[1,1], sort(c(dat[dif,1], dat[dif+1,1])), dat[nrow(dat),1])
      dif <- matrix(dif, ncol=2, byrow=TRUE)
      probs <- c()
      for(i in 1:nrow(dif))
        probs[i] <- round(100*sum(dat[which(dat[,1]==dif[i,1]):which(dat[,1]==dif[i,2]),2]), 1)
      hpds <- cbind(dif, probs)
    }
  hpds
}
