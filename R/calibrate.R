#' @title Find the calibrated distributions of 14C dates
#' 
#' @param cage numeric. The uncalibrated age.
#' @param error numeric.  The lab error on the radiocarbon date.
#' @param reservoir numeric.  Resevoir effect.
#' @param prob numeric. Value from 0 - 1, the significance level.
#' @param cc numeric.  Which calibration curve should be used [0 - 5].
#' @param cc1="IntCal13.14C"
#' @param cc2="Marine13.14C"
#' @param cc3="SHCal13.14C"
#' @param cc4="mixed.14C"
#' @param cc5="gluedHemispheres.14C"
#' @param postbomb boolean.  Do we use the post-bomb curve?
#' @param pb1="postbomb_NH1.14C"
#' @param pb2="postbomb_NH2.14C"
#' @param pb3="postbomb_NH3.14C"
#' @param pb4="postbomb_SH1-2.14C"
#' @param pb5="postbomb_SH3.14C"
#' @param yrsteps=1
#' @param pbsteps=0.01
#' @param hpdsteps=1
#' @param calibt=FALSE
#' @param yrmin=c()
#' @param yrmax=c()
#' @param minC14=c()
#' @param maxC14=c()
#' @param times=5
#' @param calheight=0.3
#' @param expand=0.1
#' @param threshold=1e-6
#' @param storedat=FALSE
#' @param graph=TRUE
#' @param xlab=c()
#' @param ylab=c()
#' @param BCAD=FALSE
#' @param mar=c(3.5,3,2,1)
#' @param mgp=c(2,1,0)
#' @param bty="l"
#' @param title=c()
#' @param date.col="red"
#' @param cc.col=rgb(0,.5,0,0.7)
#' @param dist.col=rgb(0,0,0,0.3)
#' @param sd.col=rgb(0,0,0,0.5)
#' @export

calibrate <- function(cage=2450, error=50, reservoir=0, prob=0.95, cc=1, cc1="IntCal13.14C", cc2="Marine13.14C", cc3="SHCal13.14C", cc4="mixed.14C", cc5="gluedHemispheres.14C", postbomb=FALSE, pb1="postbomb_NH1.14C", pb2="postbomb_NH2.14C", pb3="postbomb_NH3.14C", pb4="postbomb_SH1-2.14C", pb5="postbomb_SH3.14C", yrsteps=1, pbsteps=0.01, hpdsteps=1, calibt=FALSE, yrmin=c(), yrmax=c(), minC14=c(), maxC14=c(), times=5, calheight=0.3, expand=0.1, threshold=1e-6, storedat=FALSE, graph=TRUE, xlab=c(), ylab=c(), BCAD=FALSE, mar=c(3.5,3,2,1), mgp=c(2,1,0), bty="l", title=c(), date.col="red", cc.col=rgb(0,.5,0,0.7), dist.col=rgb(0,0,0,0.3), sd.col=rgb(0,0,0,0.5))
  .calibrate(cage, error, reservoir, prob, cc, cc1, cc2, cc3, cc4, cc5, postbomb, pb1, pb2, pb3, pb4, pb5, yrsteps, pbsteps, hpdsteps, calibt, yrmin, yrmax, minC14, maxC14, times, calheight, expand, threshold, storedat, graph, xlab, ylab, BCAD, mar, mgp, bty, title, date.col, cc.col, dist.col, sd.col)

.calibrate <- function(cage, error, reservoir, prob, cc, cc1, cc2, cc3, cc4, cc5, postbomb, pb1, pb2, pb3, pb4, pb5, yrsteps, pbsteps, hpdsteps, calibt, yrmin, yrmax, minC14, maxC14, times, calheight, expand, threshold, storedat, graph, xlab, ylab, BCAD, mar, mgp, bty, title, date.col, cc.col, dist.col, sd.col)
{
  # set calibration curve
  if(cc==1) calcurve <- read.table(cc1) else
    if(cc==2) calcurve <- read.table(cc2) else
      if(cc==3) calcurve <- read.table(cc3) else
        if(cc==4) calcurve <- read.table(cc4) else
          if(cc==5) calcurve <- read.table(cc5) else
            stop("I do not understand which calibration curve you mean, check the manual", call.=FALSE)
  
  # include postbomb curve if required
  if(cage < 0)
  {
    pb <- 0
    if(postbomb==FALSE)
      stop("\n  Negative 14C age, should I use a postbomb curve?\n", call.=FALSE)
    if(postbomb==1) pb <- pb1 else
      if(postbomb==2) pb <- pb2 else
        if(postbomb==3) pb <- pb3 else
          if(postbomb==4) pb <- pb4 else
            if(postbomb==5) pb <- pb5 else
              stop("I do not understand which postbomb curve you mean, check the manual", call.=FALSE)
    yrsteps <- min(pbsteps, yrsteps)
    if(length(pb) > 0)
    {
      pb <- read.table(pb)
      pb.x <- seq(min(pb[,1]), max(pb[,1]), by=yrsteps)
      pb.y <- approx(pb[,1], pb[,2], pb.x)$y
      pb.sd <- approx(pb[,1], pb[,3], pb.x)$y
      calcurve <- cbind(c(pb.x, calcurve[,1]), c(pb.y, calcurve[,2]), c(pb.sd, calcurve[,3]))
    }
    cat("  postbomb date, interpolating to every", pbsteps, "yr.")
  }
  
  # check whether date lies partly or entirely beyond the calibration curve
  if(length(reservoir) == 2) # assuming that first value is mean offset, second is error
  {
    error <- sqrt(error^2 + reservoir[2]^2)
    reservoir <- reservoir[1]
  }
  border <- 0
  if(cage-reservoir-error < min(calcurve[,2]+calcurve[,3]))
    if(cage-reservoir+error > min(calcurve[,2]-calcurve[,3]))
      border <- 1 else border <- 2
  if(cage-reservoir+error > max(calcurve[,2]-calcurve[,3]))
    if(cage-reservoir-error < max(calcurve[,2]+calcurve[,3]))
      border <- 1 else border <- 2
  if(border==1)
    cat("\nDate falls partly beyond calibration curve and will be truncated!")
  if(border==2)
    stop("\nCannot calibrate dates beyond calibration curve!\n\n")
  
  # work in BC/AD if needed, and prepare for calculations in f14C
  if(BCAD)
  {
    theta <- 1950-calcurve[,1]
    ad <- max(which(theta > 0)) # one side of the border between AD and BC
    theta <- c(theta[1:(ad-1)], theta[ad]:theta[ad+2], theta[(ad+3):length(theta)])
    mu <- approx(1950-calcurve[,1], calcurve[,2], theta)$y
    sigma <- approx(1950-calcurve[,1], calcurve[,3], theta)$y
    theta[theta <=0] <- theta[theta <=0]-1
    calcurve <- cbind(theta, mu, sigma)
  } else theta <- calcurve[,1]
  
  f.mu <- exp(-calcurve[,2]/8033)
  f.sigma <- exp(-(calcurve[,2]-calcurve[,3])/8033) - f.mu
  f.cage <- exp(-(cage-reservoir)/8033)
  f.error <- f.cage - exp(-(cage-reservoir+error)/8033)
  
  # calibrate the date and report its highest posterior density (hpd) range
  if(length(xlab)==0) xlab <- ifelse(BCAD, "cal BC/AD", "cal BP")
  calib <- .caldist(f.cage, f.error, theta, f.mu, f.sigma, yrsteps, threshold, calibt, BCAD)
  hpd <- .hpd(calib, prob, hpdsteps, yrsteps)
  colnames(hpd) <- c("yrmin", "yrmax", "prob")
  dat <- list(calib=calib, hpd=hpd)
  if(storedat) dat <<- dat
  cat("\nmin\tmax\tprob\n")
  for(i in 1:nrow(hpd))
  {
    for(j in 1:3) cat(hpd[i,j], "\t")
    cat("\n")
  }
  cat("\n")
  
  # produce a graph of the calibrated distribution (default)
  if(graph)
  {
    ifelse(BCAD,
           xrange <- 1950+c((1+expand)*(min(calib[,1])-1950), (1-expand)*(max(calib[,1])-1950)),
           xrange <- c((1+expand)*max(calib[,1]), (1-expand)*min(calib[,1])))
    if(length(yrmin) > 0) xrange[2] <- yrmin
    if(length(yrmax) > 0) xrange[1] <- yrmax
    ifelse(BCAD,
           cc <- calcurve[max(which(theta >= min(xrange))):min(which(theta <= max(xrange))),],
           cc <- calcurve[min(which(theta >= min(xrange))):max(which(theta <= max(xrange))),])
    
    # first plot the calibrated distribution, and its hpd ranges
    par(mar=mar, mgp=mgp, bty=bty, xaxs="i", xaxt="s", yaxt="n", yaxs="i", new=FALSE)
    pol <- cbind(c(calib[,1], rev(calib[,1])), c(calib[,2]/max(calib[,2]), rep(0, length=nrow(calib))))
    plot(0, type="n", xlim=xrange, ylim=c(0,1/calheight), xlab="", ylab="")
    polygon(pol, col=dist.col, border=NA)
    for(i in 1:nrow(hpd))
    {
      if(hpd[i,1]==hpd[i,2])
      {
        probs <- calib[which(calib[,1]==hpd[i,1]),]
        lines(rep(probs[1], 2), c(0, probs[2]/max(calib[,2])), col=grey(.5))
      } else
      {
        probs <- calib[max(which(calib[,1]<=hpd[i,1])):max(which(calib[,1]<=hpd[i,2])),]
        pol <- cbind(c(probs[,1], rev(probs[,1])), c(probs[,2]/max(calib[,2]), rep(0, length=nrow(probs))))
        polygon(pol, col=sd.col, border=NA)
      }
    }
    lines(calib[,1], calib[,2]/max(calib[,2]))
    abline(h=0)
    
    # now draw the 14C distribution (normal distribution, on vertical axis)
    par(new=TRUE, yaxt="s", yaxs="r", xaxt="n")
    if(length(cc)==3) cc <- cbind(cc[1], cc[2], cc[3])
    if(reservoir!=0)
      main <- substitute(cage-res %+-% er, list(cage=cage, er=error, res=reservoir)) else
        main <- substitute(cage %+-% er, list(cage=cage, er=error))
    if(length(title)>0) main <- title
    if(length(minC14)==0)
      minC14 <- min(cc[,2]-qnorm(1-(1-prob)/2)*cc[,3], cage-reservoir-qnorm(1-(1-prob)/2)*error)
    if(length(maxC14)==0)
      maxC14 <- max(cc[,2]+qnorm(1-(1-prob)/2)*cc[,3], cage-reservoir+qnorm(1-(1-prob)/2)*error)
    if(length(ylab)==0)
      ylab <- expression(paste(""^14, "C BP"))
    plot(0, type="n", xlim=xrange, ylim=c(minC14, maxC14), xlab=xlab, ylab=ylab, main=main)
    if(length(calibt) > 0)
      times <- 5*times
    yage <- (cage-reservoir-times*error):(cage-reservoir+times*error) # must not be on F14C for plot
    if(length(calibt) < 2)
      xage <- dnorm(exp(-yage/8033), f.cage, f.error) else
        xage <- (calibt[2] + ((f.cage-exp(-yage/8033))^2) / (2*(f.error^2))) ^ -(calibt[1]+0.5)
    xage.plot <- xrange[1]-((xrange[1]-xrange[2])*calheight)*xage/max(xage)
    pol <- cbind(c(xage.plot, rep(xrange[1], length(xage))), c(yage, rev(yage)))
    polygon(pol, col=dist.col, border="black")
    
    # draw the highest posterior density (hpd) range of the 14C date
    xage[which(cumsum(xage)/sum(xage) > 1 - (1-prob)/2)] <- 0
    xage[which(cumsum(xage)/sum(xage) < (1-prob)/2)] <- 0
    xage <- xrange[1]-((xrange[1]-xrange[2])*calheight)*(xage/max(xage))
    pol <- cbind(c(xage, rep(xrange[1], length=length(xage))), c(yage, rev(yage)))
    polygon(pol, col=sd.col, border=FALSE)
    
    # plot the mid and error of the 14C date
    points(xrange[1]-.01*(xrange[1]-xrange[2]), cage, pch=19, col=date.col)
    lines(rep(xrange[1]-.01*(xrange[1]-xrange[2]), 2), c(cage-error, cage+error), lwd=2, col=date.col)
    
    # now draw the calibration curve
    pol <- cbind(c(theta, rev(theta)),
                 c(calcurve[,2]-qnorm(1-(1-prob)/2)*calcurve[,3], rev(calcurve[,2]+qnorm(1-(1-prob)/2)*calcurve[,3])))
    polygon(pol, border=cc.col, col=cc.col)
  }
}
