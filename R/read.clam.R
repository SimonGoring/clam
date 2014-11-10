#' Load data from the core and begin calibration.
#' 
#' Read the data and perform first calculations incl. calibrations
#' 
#' @param name
#' @param ext
#' @param hpdsteps
#' @param yrsteps
#' @param prob
#' @param times
#' @param sep
#' @param BCAD
#' @param storedat
#' @param ignore
#' @param thickness
#' @param youngest
#' @param slump
#' @param threshold
#' @param theta
#' @param f.mu
#' @param f.sigma
#' @param calibt
#' @param extradates
#' @param calcurve
#' @param postbomb
#' 
#' @return The age model. (?)
#' @export
#' 
.read.clam <- function(name, ext, hpdsteps, yrsteps, prob, times, sep, BCAD, storedat, ignore, thickness, youngest, slump, threshold, theta, f.mu, f.sigma, calibt, extradates, calcurve, postbomb)
{
  # read the file with the dating information
  dat <- list(coredir=paste("Cores/", name, "/", sep=""), name=name)
  if(!file.exists(paste("Cores/", name, sep="")))
    stop(paste("\n\n Warning, cannot find a folder within Cores/ named ", name, ". Have you saved it in the right place and with the right name? Please check the manual\n\n", sep=""), call.=FALSE)
  if(!file.exists(paste(dat$coredir, name, ext, sep="")))
    stop(paste(" \n\n Warning, cannot find file ", name, ".csv in folder Cores/", name, ". Have you saved it in the right place and named it correctly? Please check the manual\n\n", sep=""), call.=FALSE)
  dets <- suppressWarnings(read.table(paste(dat$coredir, name, ext, sep=""), comment.char="", header=TRUE, sep=sep, na.strings = c("#N/A!", "NA", "@NA")))
  
  # ignore dates if required, add thickness column if it was left out
  if(length(ignore) > 0)
  {
    dat$ignore <- as.character(dets[ignore,1])
    dets <- dets[-ignore,]
  }
  if(ncol(dets) < 7)
    dets <- cbind(dets, thickness) else
      dets[is.na(dets[,7]),7] <- thickness
  
  # should slumps be taken into account?
  if(length(slump) > 0)
  {
    d.adapt <- dets[,6]
    d.lost <- c()
    for(i in 1:nrow(slump))
    {
      below.slump <- which(dets[,6] > max(slump[i,]))
      above.slump <- which(dets[,6] < min(slump[i,]))
      d.lost <- c(d.lost, which(!(1:nrow(dets) %in% c(above.slump, below.slump))))
      d.adapt[below.slump] <- d.adapt[below.slump] - (max(slump[i,])-min(slump[i,]))
    }
    dets[,6] <- d.adapt
    if(length(d.lost) > 0)
      dets <- dets[-d.lost,]
  } 
  # check for common errors
  dets <- dets[,1:7]
  x <- 0
  for(i in 2:7) if(is.factor(dets[,i])) x <- 1
  if(x==1) stop(paste("\n Some value fields in ", name, ".csv contain letters, please adapt", sep=""), call.=FALSE)
  if(length(dets[is.na(dets[,2]),2])+length(dets[is.na(dets[,3]),3]) != nrow(dets))
    stop(paste("\n Remove duplicate entries within the C14 and calendar fields in ", name, ".csv", sep=""), call.=FALSE)
  if(min(dets[,4]) <= 0)
    stop(paste("\n Errors of dates should be larger than zero. Please adapt ", name, ".csv", sep=""), call.=FALSE)
  dat$ID <- as.character(dets[,1])
  
  # correct for any reservoir effect
  dets[is.na(dets[,5]),5] <- 0
  dat$cage <- dets[,2] - dets[,5]
  dat$error <- dets[,4]
  
  # work in F14C for calibration
  dat$f.cage <- exp(-dat$cage/8033)
  dat$f.error <- dat$f.cage - exp(-(dat$cage+dat$error)/8033)
  
  # check if any 14C dates are (entirely or partly) beyond the calibration curve
  outside <- which(!is.na(dat$cage))
  rangecc <- c(min(calcurve[,2]-calcurve[,3]),max(calcurve[,2]+calcurve[,3]))
  outside <- outside[c(which(dat$cage[outside]-times*dat$error[outside] < rangecc[1]), which(dat$cage[outside]+times*dat$error[outside] > rangecc[2]))]
  if(length(outside) > 0)
  {
    truncate <- 0
    for(i in 1:length(outside)) # check if date lies only partly beyond the curve limits
      if((dat$cage[outside[i]]-times*dat$error[outside[i]] < rangecc[1] &&
            dat$cage[outside[i]]+times*dat$error[outside[i]] > rangecc[1]) ||
           (dat$cage[outside[i]]-times*dat$error[outside[i]] < rangecc[2] &&
              dat$cage[outside[i]]+times*dat$error[outside[i]] > rangecc[2]))
        truncate <- truncate + 1
    if(truncate > 0)
      cat("\n Warning, dates spanning beyond the calibration curve will be truncated! ")
    
    # remove dates which lie entirely outside the limits of the calibration curve
    outside <- outside[c(which(dat$cage[outside]+qnorm(1-(1-prob)/2)*dat$error[outside] < rangecc[1]), which(dat$cage[outside]-qnorm(1-(1-prob)/2)*dat$error[outside] > rangecc[2]))]
    if(length(outside) > 0)
    {
      cat("\n Warning, dates older than the calibration curve will be ignored! ")
      dets <- dets[-outside,]
      dat$cage <- dat$cage[-outside]
      dat$error <- dat$error[-outside]
      dat$f.cage <- dat$f.cage[-outside]
      dat$f.error <- dat$f.error[-outside]
      dat$outside <- dat$ID[outside]
      dat$ID <- dat$ID[-outside]
    }
  }
  
  # fill the 'dat' list with additional information
  dat$cal <- c(dets[,3], extradates)
  dat$res <- c(dets[,5], extradates)
  dat$depth <- c(dets[,6], extradates)
  dat$thick <- c(dets[,7], rep(thickness, length(extradates)))
  dat$BCAD <- BCAD
  
  # find distribution (calibrated if 14C) and point estimates for each date
  for(i in 1:length(dat$depth))
  {
    if(length(extradates) > 0 && i > nrow(dets))
    {
      tmp <- read.table(paste(dat$coredir, name, "_", extradates[i-nrow(dets)], ".txt", sep=""))
      calib <- cbind(tmp[,1], tmp[,2]/sum(tmp[,2]))
    } else
      if(is.na(dat$cage[[i]]))
      {
        age <- dat$cal[[i]]
        error <- dat$error[[i]]
        ageseq <- seq(age-(times*error), age+(times*error), by=yrsteps)
        calib <- cbind(ageseq, dnorm(ageseq, age, error))
      } else
        calib <- .caldist(dat$f.cage[[i]], dat$f.error[[i]], theta, f.mu, f.sigma, yrsteps, threshold, calibt, BCAD)
    if(length(youngest) > 0) # truncate ages younger than a limit
    {
      if(BCAD) calib <- calib[which(calib[,1] <= youngest),] else
        calib <- calib[which(calib[,1] >= youngest),]
      if(length(calib) == 0)
        if(BCAD)
          calib <- cbind(seq(youngest-(3*yrsteps), youngest+yrsteps, length=5), c(0:3,0)/3) else
            calib <- cbind(seq(youngest-yrsteps, youngest+(3*yrsteps), length=5), c(0,3:0)/3)
    }
    dat$calib[[i]] <- calib
    dat$hpd[[i]] <- .hpd(calib, prob=prob, hpdsteps=hpdsteps, yrsteps=yrsteps)
    dat$mid1[[i]] <- (dat$hpd[[i]][1] + dat$hpd[[i]][2*nrow(dat$hpd[[i]])])/2
    yrs <- calib[,1]
    dat$mid2[[i]] <- mean(c(max(yrs), min(yrs)))
    dat$wmn[[i]] <- weighted.mean(calib[,1], 1/calib[,2])
    dat$med[[i]] <- calib[max(which(cumsum(calib[,2]) <= .5)),1]
    dat$mode[[i]] <- calib[which(calib[,2] == max(calib[,2])),1][1]
  }
  
  if(storedat) dets <<- dets
  dat
}
