#' Mixed effects modelling using calibrated dates.
#' 
#' Akin to Heegaard et al.'s mixed effect modelling, but using calibrated dates.  This function is internal to the package.
#' @param its numeric.  The number of iterations used.
#' @param depths
#' @param cals
#' @param cages
#' @param errors
#' @param calibs
#' @param Est
#' @param theta
#' @param f.mu
#' @param f.sigma
#' @param yrsteps
#' @param calibt
#' 
#' @return A vector of ages.
#' @export

.mixed.effect <- function(its, depths, cals, cages, errors, calibs, Est, theta, f.mu, f.sigma, yrsteps, calibt)
{
  cat("\n Mixed effect modelling, this will take some time")
  smp <- array(1, dim=c(length(depths), 1+its, 2))
  smp[,1,1] <- Est
  for(i in 1:length(cals))
    if(!is.na(cals[i]))
      if(length(calibt)==0)
      {
        x <- rnorm(its, cals[i], errors[i])
        smp[i,(1:its)+1,] <- c(x, dnorm(x, cals[i], errors[i]))
      } else
      {
        x <- (cals[i]-10*errors[i]) : (cals[i]+10*errors[i])
        x <- cbind(x, .calibt(calibt[1], calibt[2], cals[i], errors[i], x, x, 0))
        o <- order(x[,2], decreasing=TRUE)
        x <- cbind(x[o,1], cumsum(x[o,2])/sum(x[,2]))
        sampled.x <- max(which(x[,2] <= runif(1, 0, max(x[,2]))))
        smp[i,(1:its)+1,] <- x[sampled.x,]
      } else
        for(j in 1:its)
        {
          if(j/(its/3) == round(j/(its/3))) cat(".")
          yr <- rnorm(1, cages[i], errors[i])
          f.yr <- exp(-yr/8033)
          f.error <- f.yr - exp(-(yr+errors[i])/8033)
          yr <- cbind(theta, dnorm(f.mu, f.yr, sqrt(f.error^2+f.sigma^2)))
          yr <- yr[yr[,2]>0,]
          yr <- approx(yr[,1], yr[,2], seq(min(yr[,1]), max(yr[,1]), by=yrsteps))
          smp.yr <- sample(length(yr$x), 1, prob=yr$y)
          smp[i,j+1,] <- c(yr$x[smp.yr], yr$y[smp.yr])
        }
  smp
}
