#' @title Internal function for age-model plotting.
#' 
#' @param yrmin numeric.
#' @param yrmax numeric.
#' @param dmin numeric.
#' @param dmax numeric.
#' @param revaxes boolean.
#' @param revd boolean.
#' @param revyr boolean.
#' @param yrlab
#' @param dlab
#' @param hiatus
#' @param depthseq
#' @param outliers
#' @param plotrange
#' @param BCAD,
#' @param greyscale
#' @param chron
#' @param C14col, 
#' @param outcol
#' @param outlsize
#' @param bestcol
#' @param rangecol
#' @param dat
#' @param calrange
#' @param depth
#' @param calhght
#' @param maxhght
#' @param mirror
#' @param calcol
#' @param slump
#' @param slumpcol
#' @param plotname
#' @param name
#' @param bty="l"
#' @param mar
#' @param mgp
#' @param ash=FALSE
#' 
#' @return The age model plot.
#' @export

.ageplot <- function(yrmin, yrmax, 
                     dmin, dmax, 
                     revaxes, revd, revyr, 
                     yrlab, dlab, 
                     hiatus, depthseq, outliers, plotrange, 
                     BCAD, 
                     greyscale, chron, C14col, 
                     outcol, outlsize, 
                     bestcol, rangecol, 
                     dat, calrange, 
                     depth, calhght, maxhght, 
                     mirror, calcol, slump, slumpcol, 
                     plotname, name, bty="l", mar, mgp, ash=FALSE)
{
  # set up initial parameters
  if(length(dlab)==0) dlab <- paste("Depth (", depth, ")", sep="")
  ifelse(BCAD || !revyr, yr.lim <- c(yrmin, yrmax), yr.lim <- c(yrmax, yrmin))
  if(revd) d.lim <- c(dmax, dmin) else d.lim <- c(dmin, dmax)
  
  par(xaxt="s", xaxs="r", yaxt="s", yaxs="r", bty=bty, mar=mar, mgp=mgp, font=2)
  if(revaxes){
    plot(0, type="n", ylim=yr.lim, xlim=d.lim, xlab=dlab, ylab=yrlab)
  }
  if(!revaxes){
    plot(0, type="n", xlim=yr.lim, ylim=d.lim, xlab=yrlab, ylab=dlab)
  }
  
  if(plotname) legend("topleft", name, bty="n")
  
  # draw histograms of all age-depth models. Off by default, time-consuming!
  if(length(greyscale)==1)
  {
    plotrange <- FALSE
    depgr=seq(dmin, dmax, length=greyscale)
    for(i in 2:greyscale)
    {
      temp <- density(chron[max(which(calrange[,1]<=depgr[i])),2:ncol(chron)], n=greyscale)
      if(revaxes){
        image(c(depgr[i-1], depgr[i]), temp$x, matrix(temp$y), col=grey(1-(0:100)/100), add=TRUE)
      }
      if(!revaxes){
        image(temp$x, c(depgr[i-1], depgr[i]), matrix(temp$y), col=grey(1-(0:100)/100), add=TRUE)
      }
    }
  }
  
  # draw the age-depth models, per section if hiatuses were inferred
  if(length(hiatus) > 0)
  {
    if(length(slump) == 0)
      hiatusseq <- sort(c(range(depthseq), hiatus)) else
        hiatusseq <- sort(c(range(depthseq, depthseq+sum(slump[,2]-slump[,1])), hiatus))
    for(i in 2:length(hiatusseq))
    {
      sec <- calrange[min(which(calrange[,1] > hiatusseq[i-1])):max(which(calrange[,1] < hiatusseq[i])),]
      pol <- cbind(c(sec[,2], rev(sec[,3])), c(sec[,1], rev(sec[,1])))
      if(plotrange)
        if(revaxes)
          polygon(pol[,2], pol[,1], col=rangecol, border=rangecol) else
            polygon(pol, col=rangecol, border=rangecol)
      if(revaxes)
        lines(sec[,1], sec[,4], lwd=2, col=bestcol) else
          lines(sec[,4], sec[,1], lwd=2, col=bestcol)
      if(revaxes)
        abline(v=hiatus, col="grey", lty="dashed") else
          abline(h=hiatus, col="grey", lty="dashed")
    }
  } else
  {
    pol <- cbind(c(calrange[,2], rev(calrange[,3])), c(calrange[,1], rev(calrange[,1])))
    if(plotrange)
      if(revaxes)
        polygon(pol[,2], pol[,1], col=rangecol, border=rangecol) else
          polygon(pol, col=rangecol, border=rangecol)
    if(revaxes)
      lines(calrange[,1], calrange[,4], lwd=2, col=bestcol) else
        lines(calrange[,4], calrange[,1], lwd=2, col=bestcol)
  }
  
  # draw slumps if these were given
  if(length(slump) > 0)
    for(i in 1:nrow(slump))
      if(revaxes)
        rect(min(slump[i,]), min(yr.lim)-1e4, max(slump[i,]), max(yr.lim)+1e4, col=slumpcol, border=slumpcol) else
          rect(min(yr.lim)-1e4, min(slump[i,]), max(yr.lim)+1e4, max(slump[i,]), col=slumpcol, border=slumpcol)
  
  # draw the calibrated distributions of the dates
  top <- 1
  for(i in 1:length(dat$depth))
    top <- min(top, max(dat$calib[[i]][,2])) # find the lowest peak
  
  if(calhght > 0)
    for(i in 1:length(dat$depth))
    {
      if(is.na(dat$cal[[i]])) col <- C14col else col <- calcol
      pol <- dat$calib[[i]] # already normalised to 1
      if(ash) pol[,2] <- pol[,2]/max(pol[,2])/1e3 # draw all same height
      pol[pol[,2] > maxhght,2] <- maxhght
      pol[,2] <- calhght*(dmax-dmin)*pol[,2]/(top*100)
      pol <- cbind(c(pol[,1], rev(pol[,1])),
                   c(dat$depth[[i]]-pol[,2], dat$depth[[i]]+mirror*rev(pol[,2])))
      if(revaxes) polygon(pol[,2], pol[,1], col=col, border=col) else
        polygon(pol, col=col, border=col)
    }
  
  # draw the calibrated ranges of the dates
  for(i in 1:length(dat$depth))
  {
    if(is.na(dat$cal[[i]])) col <- C14col else col <- calcol
    for(j in 1:nrow(dat$hpd[[i]]))
      if(revaxes)
        rect(dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,1], dat$depth[i]+dat$thick[i]/2, dat$hpd[[i]][j,2], lwd=1, lend=2, col=col, border=NA) else
          rect(dat$hpd[[i]][j,1], dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,2], dat$depth[i]+dat$thick[i]/2, lwd=1, lend=2, col=col, border=NA)
  }
  if(length(outliers)>0) # any outliers?
  {
    for(i in outliers)
      for(j in 1:nrow(dat$hpd[[i]]))
        if(revaxes)
          rect(dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,1], dat$depth[i]+dat$thick[i]/2, dat$hpd[[i]][j,2], col=outcol, border=outcol, lwd=1, lend=2) else
            rect(dat$hpd[[i]][j,1], dat$depth[i]-dat$thick[i]/2, dat$hpd[[i]][j,2], dat$depth[i]+dat$thick[i]/2, col=outcol, border=outcol, lwd=1, lend=2)
    if(revaxes)
      points(dat$depth[outliers], dat$mid1[outliers], cex=outlsize, pch=4, col=outcol) else
        points(dat$mid1[outliers], dat$depth[outliers], cex=outlsize, pch=4, col=outcol)
  }
}
