#' write files of the age-depth model, calibrated ranges, and settings
#' @param dat
#' @param runname
#' @param calrange
#' @param name
#' @param prob
#' @param type
#' @param remove.reverse
#' @param smooth
#' @param wghts
#' @param its
#' @param outliers
#' @param ignore
#' @param est
#' @param BCAD
#' @param yrsteps
#' @param every
#' @param decimals
#' @param accrate
#' @param depth
#' @param depthseq
#' @param hiatus
#' @param gfit
#' @param reversal
#' @param plotpdf
#' @param plotpng
#' @param yrmin
#' @param yrmax
#' @param dmin
#' @param dmax
#' @param yrlab
#' @param dlab
#' @param plotrange
#' @param greyscale
#' @param chron
#' @param C14col
#' @param outcol
#' @param outlsize
#' @param bestcol
#' @param rangecol
#' @param calhght
#' @param maxhght
#' @param mirror
#' @param calcol
#' @param slump
#' @param slumpcol
#' @param revaxes
#' @param revyr
#' @param revd
#' @param calibt
#' @param youngest
#' @param extradates
#' @param plotname
#' @param calcurve
#' @param ccname
#' @param postbomb
#' @param pbnames
#' @param depths.file
#' @param bty
#' @param mar
#' @param mgp
#' @param ash 
#' 
#' @return Something
#' @export
#' 
.write.clam <- function(dat, runname, calrange, name, prob, type, remove.reverse, smooth, wghts, its, outliers, ignore, est, BCAD, yrsteps, every, decimals, accrate, depth, depthseq, hiatus, gfit, reversal, plotpdf, plotpng, yrmin, yrmax, dmin, dmax, yrlab, dlab, plotrange, greyscale, chron, C14col, outcol, outlsize, bestcol, rangecol, calhght, maxhght, mirror, calcol, slump, slumpcol, revaxes, revyr, revd, calibt, youngest, extradates, plotname, calcurve, ccname, postbomb, pbnames, depths.file, bty, mar, mgp, ash)
{
  # age-depth model; age estimates, accumulation rates and ranges for every analysed depth
  runnames <- c("_interpolated", "_polyn_regr", "_cubic_spline", "_smooth_spline", "_loess")

  calrange <- cbind(calrange, round(c(diff(calrange[,4])/diff(calrange[,1]), NA), decimals+2))
  
  if(accrate==1) calrange[,5] <- 1/calrange[,5]
  
  calrange[,2:4] <- round(calrange[,2:4], decimals)
  
  ifelse(length(runname)==0, runname <- runnames[type], runname)
  
  if(depths.file && file.exists(dd <- paste("Cores/", name, "/", name, "_depths.txt", sep="")))
  {
    dd <- read.table(dd)[,1]
    this <- c()
    for(i in 1:length(dd))
      this[i] <- which(calrange[,1]==dd[i])[1] # find where the relevant ages are
    write.table(calrange[this,], paste(dat$coredir, name, runname, "_ages.txt", sep=""), row.names=FALSE, col.names=c("depth", paste("min", 100*prob, "%", sep=""), paste("max", 100*prob, "%", sep=""), "best", "acc.rate"), quote=FALSE, sep="\t")
  } else
    write.table(calrange, paste(dat$coredir, name, runname, "_ages.txt", sep=""), row.names=FALSE, col.names=c("depth", paste("min", 100*prob, "%", sep=""), paste("max", 100*prob, "%", sep=""), "best", "accrate"), quote=FALSE, sep="\t")
  
  # calibrated ranges of all dates
  hpd.file <- file(paste(dat$coredir, name, "_calibrated.txt", sep=""), "w")
  cat(paste("Calibrated age ranges at ", 100*prob, "% confidence intervals\n", sep=""), file=hpd.file)
  for(i in 1:length(dat$depth))
  {
    cat(paste("\n\nDepth: ", dat$depth[[i]], "\nyrmin\tyrmax\tprobability\n"), file=hpd.file)
    hpds <- dat$hpd[[i]]
    for(j in 1:nrow(hpds))
    {
      for(k in 1:3) cat(hpds[j,k], "\t", file=hpd.file)
      cat("\n", file=hpd.file)
    }
  }
  
  close(hpd.file)
  
  # relevant settings and results
  set.file <- file(paste(dat$coredir, name, runnames[type], "_settings.txt", sep=""), "w")
  
  cat(paste("Settings (square brackets give names of the constants)\n\n",
            "Calibration curve: ", ccname,
            if(postbomb!=FALSE) paste(",", pbnames[postbomb], "for postbomb dates"),
            "\nAge-depth model: ",
            if(type==1) "linear interpolation between dated levels [type=1]" else
              if(type==2) ifelse(length(smooth)==0, "linear regression [type=2, smooth=c()]",
                                 paste("polynomial regression [type=2] of order", smooth, "[smooth]")) else
                                   if(type==3) "cubic spline [type=3]" else
                                     if(type==4) paste("smooth spline [type=4] with spar =", ifelse(length(smooth)<1, 0.3, smooth), "[smooth]") else
                                       if(type==5) paste("locally weighted spline [type=5] with span =", ifelse(length(smooth)<1, 0.75, smooth), "[smooth]"),
            if(wghts==1) "\nWeighted by the calibrated probabilities [wghts=1]",
            if(wghts==2) "\nWeighted by the errors (1/sdev^2) [wghts=2]",
            "\nCalculations at ", 100*prob, "% confidence ranges [prob=", prob, "]",
            "\nAmount of iterations: ", its, " [its]",
            "\nCalendar age point estimates for depths based on ",
            if(est==1) "weighted average of all age-depth curves [est=1]" else
              if(est==2) "midpoints of the hpd ranges of the age-depth curves [est=2]" else
                if(est==3) "midpoints of the hpd ranges of the dated levels [est=3]" else
                  if(est==4) "weighted means of the dated levels [est=4]" else
                    if(est==5) "medians of the dated levels [est=5]" else
                      if(est==6) "modes/maxima/intercepts of the dated levels [est=6]",
            "\nCalendar scale used: ", if(BCAD) "cal BC/AD" else "cal BP",
            " [BCAD=", BCAD, "] at a resolution of ", yrsteps, " yr [yrsteps]",
            "\nAges were calculated every ", every, " [every] ", depth,
            " [depth], from ", min(depthseq), " [dmin] to ", max(depthseq), " [dmax] ", depth, sep=""), file=set.file)
  if(length(youngest) > 0) cat("\n\nDates with ages younger than", youngest, ifelse(BCAD, "BC/AD", "cal BP"), "were truncated", file=set.file)
  if(length(calibt)> 1) cat("\n\nInstead of assuming the standard Gaussian model, a student t distribution was used with t.a =", calibt[1], "and t.b =", calibt[2], "(see Christen and Perez 2009, Radiocarbon 51:1047-1059)", file=set.file)
  if(length(slump) == 2) cat("\n\nA slump was excised between", max(slump), "and", min(slump), depth, file=set.file)
  if(length(slump) > 2)
  {
    cat("\n\nSlumps were excised from ", file=set.file)
    sl <- array(sort(slump), dim=c(2, length(slump)/2))
    for(i in 1:ncol(sl))
      cat(sl[1,i], "to", sl[2,i], depth, if(i<ncol(sl)) "and ", file=set.file)
  }
  if(length(outliers) > 0)
  {
    cat("\n\nDates assumed outlying [outliers]: ", file=set.file)
    for(i in outliers) cat(i, " (", dat$ID[i], ") ", sep="", file=set.file)
  }
  if(length(ignore) > 0)
  {
    cat("\n\nDates ignored [ignore]: ", file=set.file)
    for(i in 1:length(ignore)) cat(ignore[i], " (", dat$ignore[i], ") ", sep="", file=set.file)
  }
  if(length(dat$outside) > 0)
  {
    cat("\n\nDates outside calibration curve and ignored: ", file=set.file)
    for(i in 1:length(dat$outside)) cat(dat$outside[i], " ", sep="", file=set.file)
  }
  cat(paste(
    if(length(hiatus) > 0)
      paste("\nA hiatus was inferred at", hiatus, depth, "[hiatus]"),
    "\n\nGoodness-of-fit (-log, lower is better): ", gfit,
    if(reversal) "\nSome age-depth reversals occurred"),
      if(remove.reverse) "\nAny models with age-depth reversals were removed",
      "\n\nProduced ", date(), sep="", file=set.file)
  close(set.file)
  
  if(plotpdf)
  {
    pdf(file=paste(dat$coredir, name, runname, ".pdf", sep=""))
    .ageplot(yrmin, yrmax, dmin, dmax, revaxes, revd, revyr, dlab, yrlab, hiatus, depthseq, outliers, plotrange, BCAD, greyscale, if(length(greyscale)>0) chron else c(), C14col, outcol, outlsize, bestcol, rangecol, dat, calrange, depth, calhght, maxhght, mirror, calcol, slump, slumpcol, plotname, name, bty, mar, mgp, ash)
    dev.off()
  }
  if(plotpng)
  {
    png(file=paste(dat$coredir, name, runname, ".png", sep=""))
    .ageplot(yrmin, yrmax, dmin, dmax, revaxes, revd, revyr, dlab, yrlab, hiatus, depthseq, outliers, plotrange, BCAD, greyscale, if(length(greyscale)>0) chron else c(), C14col, outcol, outlsize, bestcol, rangecol, dat, calrange, depth, calhght, maxhght, mirror, calcol, slump, slumpcol, plotname, name, bty, mar, mgp, ash)
    dev.off()
  }
  
}
