#' @title Find the calibrated distributions of 14C dates
#' 
#' @param f.cage
#' @param f.error
#' @param theta
#' @param f.mu
#' @param f.sigma
#' @param yrsteps
#' @param threshold
#' @param calibt
#' @param BCAD boolean.  Are dates expected to be in BC/AD time?
#' @param normalise boolean.  Defaults to FALSE.
#' @export

.caldist <- function(f.cage, f.error, theta, f.mu, f.sigma, yrsteps, threshold, calibt, BCAD, normalise=FALSE)
{
  if(f.cage > 1)
  {
    if(f.cage > 1) yrsteps <- min(yrsteps, .1)
    pb <- theta[which(f.mu > 1)]
    if(length(pb)==0)
      stop("help, something exploded with a postbomb date")
    x <- approx(theta, f.mu, seq(min(pb), max(pb), by=yrsteps))
    xsd <- approx(theta, f.sigma, x$x)$y
    theta <- c(x$x, theta[which(f.mu <= 0)])
    f.mu <- c(x$y, f.mu[which(f.mu <= 0)])
    f.sigma <- c(xsd, f.sigma[which(f.mu <= 0)])
    threshold <- 0
  }
  
  # calibrate; find how far f.cage (measurement) is from f.mu (calibration curve)
  if(length(calibt) < 2)
    cal <- cbind(theta, dnorm(f.mu, f.cage, sqrt(f.error^2+f.sigma^2))) else
      cal <- cbind(theta, .calibt(calibt[1], calibt[2], f.cage, f.error, theta, f.mu, f.sigma))
  
  # interpolate and normalise calibrated distribution to 1
  cal <- cal[min(which(cal[,2] > 0)):max(which(cal[,2] > 0)),] # remove unnecessary data
  cal <- approx(cal[,1], cal[,2], seq(min(cal[,1]), max(cal[,1]), by=yrsteps))
  cal <- cbind(cal$x, cal$y/sum(cal$y))
  if(BCAD && (0 %in% cal[,1]))
    cal <- cal[-which(cal[,1]==0),] # 0 BC/AD does not exist
  # only report those normalised calibrated probabilities beyond a threshold
  cal[cal[,2] > threshold,]
}
