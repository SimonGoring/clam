#' Goodness-of-fit estimates from clam models.
#' 
#' Internal function for calculating goodness of fit.
#' 
#' @param theta
#' @param f.mu
#' @param f.sigma
#' @param dat
#' @param calrange
#' @param outliers
#' 
#' @return Goodness of fit.
#' @export

.gfit <- function(theta, f.mu, f.sigma, dat, calrange, outliers)
{
  
  if(length(outliers) > 0)
  {
    dat$cage <- dat$cage[-outliers]
    dat$error <- dat$error[-outliers]
    dat$cal <- dat$cal[-outliers]
    dat$model <- dat$model[-outliers]
  }
  
  gfit <- pnorm(dat$cal, dat$model, dat$error^2)
  
  if(length(c14 <- which(!is.na(dat$cage))) > 0) # if there are radiocarbon dates
  {
    gfit.c <- approx(theta, f.mu, dat$model[c14])$y # C14 age at cc of modelled cal date
    f.cage <- exp(-dat$cage[c14]/8033)
    f.error <- exp(-(dat$cage[c14]-dat$error[c14])/8033) - f.cage
    gfit.var <- f.error^2 + approx(theta, f.sigma, dat$model[c14])$y^2
    gfit[c14] <- pnorm(f.cage, gfit.c, sqrt(gfit.var)) # deviation between measured and cc ages
  }
  
  dat$gfit <- -sum(log(gfit[!is.na(gfit)]))
}
