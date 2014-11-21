#' Age model construction.
#'
#' clam, R code for classical (non-Bayesian) age-depth modelling, this version 2.2
#' see accompanying manual.html and Blaauw 2010 (Quaternary Geochronology 5: 512-518)
#' This is the main age-depth modelling function.
#' The default values can be changed permanently within this file or temporarily when calling clam()
#' 
#' @param name character.  The site name, for use with file-based reconstructions.
#' @param type numeric.  Model type to generate the chronology.  Default is 1 (linear), also available are 2, polynomial; 3, cubic spline; 4, smooth spline; 5, locally weighted spline (loess).
#' @param smooth numeric.  Associated with the \code{type} variable.  See decription for more details.
#' @param prob numeric. Probability level for the reported confidence intervals [0 - 1].
#' @param its numeric. Number of iterations used to build the chronology.
#' @param wghts numeric. Age weighting, either one of: 1, probability weighting; 2, laboratory error.
#' @param cc numeric. Define the calibration curve to be used.  One of 1 - 5.
#' @param cc.name character. A vector of calibration curve names, 1, cc2, cc3, cc4, cc5 character.  Calibration curve data.
#' @param postbomb boolean.  Are some of the radiocarbon dates post-bomb?  Defaults to FALSE.
#' @param pb.name character.  Post-bomb calibration curve names, defaults to \code{Post + cc.name}.
#' @param outliers A single value or vector index of dates considered to be outliers (counted from the top). 
#' @param ignore A single value or vector index of dates that should be ignored entirely.
#' @param youngest A single value indicating the minimum possible date for the core.
#' @param extradates The numeric depth(s) of any extra dating information, associated with a \code{MyCore_XX.txt} file where XX is equal to \code{extradates}
#' @param slump Depth intervals to be excized from the core, listed in order from top to bottom.  Multiple slumps listed in sequence: \code{c(top1, bottom1, top2, bottom2)}
#' @param est The method used to convert the estimate to a single point value (from 1 - 7).
#' @param calibt The parameters for the t-distribution used for calibrating 14C dates, defaults to Gaussian if not provided.
#' @param mixed.effect Model error using Heegaard et al (2005) mixed effects method. Defaults to FALSE.  If TRUE the model uses a default \code{its} of 1000.
#' @param dmin Minimum depth of the model (either truncate or force extrapolation).
#' @param dmax Maximum depth of the model, as above.
#' @param every Depth resolution in cm.
#' @param yrmin Set the age-axis minimum limit for plotting (\code{xmin}).
#' @param yrmax  Set the age-axis maximum limit for plotting (\code{xmax}).
#' @param yrsteps Temporal resolution for calibration curves in years.  Default is 1.
#' @param pbsteps Temporal resolution for post-bomb calibration curves in years. Default is 0.01.
#' @param hpdsteps Temporal resolution for the highest posterior ranges in Calendar years.  Default is 1.
#' @param BCAD Are the years expressed in BC/AD (\code{TRUE}) or cal years BP (\code{FALSE}).  Default is \code{FALSE}.
#' @param decimals Number of decimals for rounding.  Default is 0, may be positive or negative.
#' @param accrate Report deposition rate (yrs/cm: \code{0}) or accumulation rate (cm/yr: \code{1})
#' @param ageofdepth Numeric or vector of depths for which ages are to be returned in message output.
#' @param depth Depth units. Default is \code{'cm'}.
#' @param depthseq Depth intervals for estimation.  Defaults to \code{seq(dmin, dmax, by = every)}.  If \code{depths.file = TRUE} will use depths file.
#' @param depths.file  If the core has an associated depths file \code{TRUE}.  Default is \code{FALSE}.
#' @param thickness  Thickness of the dated samples.
#' @param hiatus  The depths of any hiatuses.  A section must have more than two dated samples for extrapolation.
#' @param remove.reverse  Proportional threshold for model reversals (0-1) before a warning is thrown.
#' @param times The breadth of the sampling interval in standard  deviations.  Default is 5 (99.999996 pct. of CI).
#' @param sep Separator for the clam age file.  Default is \code{","}.
#' @param ext File extension for clam age file.  Default is \code{".csv"}.
#' @param runname  Used to generate unique names for individual runs.  If empty, the default is \code{"calibrate"}.
#' @param storedat  Should data be retained in the Global Environment?  Default is \code{FALSE}.  
#' @param threshold  Threshold probability for sampling exclusion fron calibration curves (default \code{1e-6})
#' @param proxies  If proxies are included (if a \code{name_proxies.csv} file exists), add them to the plot. Default \code{FALSE}.
#' @param revaxes  Reverse axes, plotting age on the x-axis and depth on the y-axis.  Default \code{FALSE}.
#' @param revd  Reverse the depth axis, plotting with 0 at the origin.  Default \code{FALSE}.
#' @param revyr  Reverse the age axis, with 0 at the origin.  Default \code{FALSE}.
#' @param calhght  Unit height for the calibrated 14C distributions. Default \code{1}.
#' @param maxhght  Maximum height for the calibrated age distributions. Default \code{0.01}.
#' @param mirror  Plot calibrated ages as violin plots or as one sided distributions. Default \code{TRUE}.
#' @param plotrange  Should the confidence ranges of the model be plotted?  Default \code{TRUE}.
#' @param bty Should the plot area be bounded? As in \code{par}.
#' @param mar Plot margins.  As in \code{par}.
#' @param mgp Axis text margins.  As in \code{par}.
#' @param plotpdf Output the plot to PDF.
#' @param plotpng Output the plot to PNG.
#' @param greyscale  Cancel the confidence interval plotting, numeric values 
#' @param yrlab Axis label for the age axis.  Default is \code{"cal BP"} or \code{"AD/BC"}.
#' @param dlab Axis label for the depth axis.  Default is \code{"depth (cm)"}.
#' @param calcol \code{rgb} fill color for the calibrated dates (polygon fill).
#' @param C14col \code{rgb} color for calibrated ranges of dates (polygon outline).
#' @param outcol Color for outlier dates.
#' @param outlsize Size of outlier dates.
#' @param bestcol Plotting color for the "best" model.
#' @param rangecol Plotting color for confidence intervals, as \code{rgb}.
#' @param slumpcol Plotting color of any slumps, as \code{rgb}.
#' @param plotname Should the core name be printed in the plot?
#' @param ash Should all distributions be plotted at the same height?
#' 
#' @author Maarten Blaauw, Simon Goring
#' @return A \code{data.frame} with depth, point and 95%CI values.  Also returns a text output and saved file to a 'Cores' folder in the working directory.
#' #' @section Note:
#' Currently the method requires you to have a \code{Cores} folder in your working directory.  This will soon be changes.
#'
#' @examples \dontrun{
#' # You need the neotoma package for this example.  If you do not have it installed, uncomment the following lines:
#' # require(devtools)
#' # install_github('neotoma', 'rOpenSci')
#' 
#' require(neotoma)
#' marion.site <- get_site(sitename='Marion Lake%')
#' marion.data <- get_dataset(siteid=marion.site$siteid, datasettype = 'pollen')
#' 
#' marion.download <- get_download(marion.data)
#' 
#' #  You need to have a 'Cores' directory in the current working directory.
#' 
#' if(!'Cores' %in% list.files(include.dirs=TRUE)){
#'   dir.create('Cores')
#' }
#' 
#' write_agefile(marion.download[[1]],chronology = 1,path = '.', corename = 'Marion', cal.prog = 'Clam')
#' 
#' #  Build a model using a smooth spline:
#' marion.ages <- clam('Marion', type = 4, 
#'                     depthseq=marion.download[[1]]$sample.meta$depth)
#' 
#' # plot the difference between the Neotoma model and the new clam model:
#' plot(marion.download[[1]]$chronologies$'NAPD 1'[,2], marion.ages[,4],
#'      xlab = 'NAPD Ages (uncalibrated)', ylab = 'Clam Ages')
#' abline(0,1)
#' }
#' @references
#'  Blaauw, M., 2010. Methods and code for 'classical' age-modelling of radiocarbon sequences. Quaternary Geochronology 5, 512-518
#' @keywords IO connection
#' @export

clam <- function(name="Example", type=1, smooth=c(), prob=0.95, its=1000, wghts=1, 
                 cc=1, cc.name = c("IntCal13.14C", "Marine13.14C", "SHCal13.14C", "mixed.14C", "gluedHemispheres.14C"), 
                 postbomb=FALSE, pb.name = c("postbomb_NH1.14C", "postbomb_NH2.14C", "postbomb_NH3.14C", "postbomb_SH1-2.14C", "postbomb_SH3.14C"), 
                 outliers=c(), ignore=c(), youngest=c(), extradates=c(), slump=c(), 
                 est=1, calibt=FALSE, mixed.effect=FALSE, dmin=c(), dmax=c(), every=1, 
                 yrmin=c(), yrmax=c(), yrsteps=1, pbsteps=0.01, hpdsteps=1, BCAD=FALSE, 
                 decimals=0, accrate=0, ageofdepth=c(), depth="cm", depthseq=c(), 
                 depths.file=FALSE, thickness=1, hiatus=c(), remove.reverse=0.5, 
                 times=5, sep=",", ext=".csv", runname=c(), storedat=FALSE, 
                 threshold=1e-6, proxies=FALSE, revaxes=FALSE, revd=TRUE, revyr=TRUE, 
                 calhght=0.3, maxhght=0.01, mirror=TRUE, plotrange=TRUE, bty="l", 
                 mar=c(3.5,3,2,1), mgp=c(2,1,0), plotpdf=TRUE, plotpng=TRUE, 
                 greyscale=c(), yrlab=c(), dlab=c(), calcol=rgb(0,0.5,0.5,0.5), 
                 C14col=rgb(0,0,1,0.5), outcol="red", outlsize=1, bestcol="black", 
                 rangecol=rgb(0,0,0,0.3), slumpcol=grey(0.75), plotname=TRUE, ash=FALSE)

  .clam(name, type, smooth, prob, its, wghts, cc, cc.name, postbomb, pb.name, outliers, ignore, youngest, extradates, slump, est, calibt, mixed.effect, dmin, dmax, every, yrmin, yrmax, yrsteps, pbsteps, hpdsteps, BCAD, decimals, accrate, ageofdepth, depth, depthseq, depths.file, thickness, hiatus, remove.reverse, times, sep, ext, runname, storedat, threshold, proxies, revaxes, revd, revyr, calhght, maxhght, mirror, plotrange, bty, mar, mgp, plotpdf, plotpng, greyscale, yrlab, dlab, 
        calcol, C14col, outcol, outlsize, bestcol, rangecol, slumpcol, plotname, ash)


#' @export
.clam <- function(name, type, smooth, prob, its, wghts, cc, cc.name, postbomb, pb.name, 
                  outliers, ignore, youngest, extradates, slump, est, calibt, 
                  mixed.effect, dmin, dmax, every, yrmin, yrmax, yrsteps, 
                  pbsteps, hpdsteps, BCAD, decimals, accrate, ageofdepth, depth, 
                  depthseq, depths.file, thickness, hiatus, remove.reverse, 
                  times, sep, ext, runname, storedat, threshold, proxies, 
                  revaxes, revd, revyr, calhght, maxhght, mirror, plotrange, bty, 
                  mar, mgp, plotpdf, plotpng, greyscale, yrlab, dlab, calcol, C14col, 
                  outcol, outlsize, bestcol, rangecol, slumpcol, plotname, ash) {
  
  is.between <- function(x, a, b) {
    # A small function to clean up all the || arguments below.
    (x - a)  *  (b - x) >= 0
  }

# Are there any abnormal settings provided?
  abnormal <- c(!is.between(type, 1, 6), 
                !is.between(prob, 0, 1), 
                its < 100,
                !is.between(wghts, 0, 1),
                !is.between(est, 1, 7),
                yrsteps <= 0,
                hpdsteps <= 0,
                every <= 0,
                decimals < 0,
                !is.between(cc, 1, 5),
                !is.between(accrate, 0, 1),
                thickness < 0,
                times < 1,
                calhght < 0,
                (type==5 && length(hiatus)>0))
  if(any(abnormal)){
    cat(abnormal)
    stop("\n Warning, clam cannot run with these settings! Please check the manual.\n\n", call.=FALSE)
  }

  if(name == 'Example'){
    if('Cores/Example' %in% list.files(recursive=TRUE)){
      name <- 'Example'
    }
    else{
      dir.create('Cores/Example', showWarnings = FALSE, recursive = TRUE)
      file.example <- system.file("extdata", "Example/Example.csv", package="clam")
      
      file.copy(file.example, 'Cores/Example')
    }
  }
  
  hold.out <- paste('Cores/', name, sep = '')
  cores.exists <- any(regexpr(hold.out, list.files(full.names=TRUE, recursive=TRUE), perl=TRUE) > 0)
  
  if(!cores.exists){
    stop('\nclam requires the presence of a Cores directory.\nMake sure your working directory is set properly (e.g., using setwd), and that a file exists for your core.')
  }
  dets <- suppressWarnings(read.csv(paste("Cores/", name, "/", name, ext, sep=""), sep=sep))
  d <- dets[,6]
  
  if(min(diff(d)) < 0) {
    cat("\n Warning, depths not in ascending order (top ones should come first).\n\n")
  }

  # avoid confusing warning when using sample for first time in session
  tmp <- suppressWarnings(sample(1:1e3, 1, prob=rep(.001,1e3), replace=TRUE))

  # avoid Windows/Mac habit of silently adding .txt extension to plain text files
  win <- list.files(paste("Cores/", name, sep=""), pattern=".csv.txt")
  
  if(length(win) > 0)
  {
    cat("\nRemoving unnecessary .txt extension from .csv file", win[1], "\n")
    file.rename(paste("Cores/", name, "/", name, ".csv.txt", sep=""),
    paste("Cores/", name, "/", name, ".csv", sep=""))
  }
  
  # set the calibration curve
  calcurve <- read.table(system.file("extdata", cc.name[cc], package="clam"))
  
  ccname <- cc.name[cc]
  
  # negative C14 ages and postbomb curve
  # pb <- 0
  pbnames <- pb.name
  cdat <- dets[,2]
  if(length(cdat[!is.na(cdat)]) > 0)
    if(min(cdat[!is.na(cdat)]) < 0)
      if(postbomb==FALSE)
        cat("\n\n  Warning, negative 14C ages, should I use a postbomb curve?\n") else
          {
            if(postbomb>5)
              stop("I do not understand which postbomb curve you mean, check the manual", call.=FALSE)
            yrsteps <- min(pbsteps, yrsteps)
            pb <- read.table(system.file("extdata", pbnames[postbomb], package="clam"))
            pb.x <- seq(min(pb[,1]), max(pb[,1]), by=yrsteps)
            pb.y <- approx(pb[,1], pb[,2], pb.x)$y
            pb.sd <- approx(pb[,1], pb[,3], pb.x)$y
            calcurve <- cbind(c(pb.x, calcurve[,1]), c(pb.y, calcurve[,2]), c(pb.sd, calcurve[,3]))
          }

  # work in BC/AD if needed, and prepare for calculations in f14C
  if(BCAD)
    {
      theta <- 1950-calcurve[,1]
      border <- max(which(theta > 0))
      theta <- c(theta[1:(border-1)], theta[border]:theta[border+2], theta[(border+3):length(theta)])
      mu <- approx(1950-calcurve[,1], calcurve[,2], theta)$y
      sigma <- approx(1950-calcurve[,1], calcurve[,3], theta)$y
      theta[theta <=0] <- theta[theta <=0]-1
      calcurve <- cbind(theta, mu, sigma)
    } else theta <- calcurve[,1]
  if(length(yrlab)==0) yrlab <- ifelse(BCAD, "cal BC/AD", "cal BP")
  f.mu <- exp(-calcurve[,2]/8033)
  f.sigma <- exp(-(calcurve[,2]-calcurve[,3])/8033) - f.mu

  # prepare for slumps and hiatuses
  cat(paste("Core name:", name))
  if(length(greyscale) > 0) storedat <- TRUE
  if(length(slump) > 0) {
      if(length(slump) %% 2 == 1) {
        stop("\n Warning, slumps need both upper and lower depths. Please check the manual", call.=FALSE)
      }
      
      slump <- matrix(sort(slump), ncol=2, byrow=TRUE)
      
      if(length(dmax)==0) {
        dmax <- max(dets[,6])
      }
      
      if(length(extradates) > 0) {
        dmax <- max(dmax, extradates)
      }
      
      for(i in 1:nrow(slump)) {
          d[d > min(slump[i,])] <- d[d > min(slump[i,])] - (max(slump[i,]) - min(slump[i,]))
          dmax <- dmax - (max(slump[i,])-min(slump[i,]))
      }
      if(length(hiatus) > 0)
        for(i in 1:nrow(slump))
          {
            below.slump <- which(hiatus > max(slump[i,]))
            above.slump <- which(hiatus < min(slump[i,]))
            hiatus[below.slump] <- hiatus[below.slump] - (max(slump[i,])-min(slump[i,]))
            hiatus <- hiatus[c(above.slump, below.slump)]
          }
    }

  # read in the data
  dat <- .read.clam(name, ext, hpdsteps, yrsteps, prob, times, sep, BCAD, storedat, ignore, thickness, youngest, slump, threshold, theta, f.mu, f.sigma, calibt, extradates, calcurve, pb)
  cat("\n Calibrating dates... ")
  
  # calculate the depths to be used, based on the ranges and resolution
  if(length(dmin)==0){
    dmin <- floor(min(dat$depth))
  }
  if(length(dmax)==0){
    dmax <- ceiling(max(dat$depth))
  }
  
  if(depths.file){
    #  If there is a set of depths at which we are interested in ages:
    if(file.exists(dd <- paste("Cores/", name, "/", name, "_depths.txt", sep=""))){
      if(length(depthseq) == 0){
        depthseq <- seq(dmin, dmax, by=every)
      }
      
      depthseq <- sort(unique(c(depthseq, suppressWarnings(read.table(dd))[,1])))
      dmin <- min(depthseq)#, read.table(dd)[,1])
      dmax <- max(depthseq)#, read.table(dd)[,1])
    } 
    else{
      stop(paste("\nCannot find file ", dat$name, "_depths.txt!\n", sep=""), call.=FALSE)
    }
  }
  
  if(length(depthseq) == 0){
    depthseq <- seq(dmin, dmax, by=every)
  }
  
  if(proxies) {
    storedat <- TRUE
    
    if(file.exists(dd <- paste("Cores/", name, "/", name, "_proxies.csv", sep=""))){
      dat$proxies <- suppressWarnings(read.csv(dd, sep=sep))
    }
    else{
      stop(paste("\nCannot find file ", dat$name, "_proxies.csv!\n", sep=""), call.=FALSE)
    }
    
    dmin <- min(depthseq, dat$proxies[,1])
    dmax <- max(depthseq, dat$proxies[,1])
    depthseq <- sort(unique(c(depthseq, dat$proxies[,1])))
  } 
      
  if(length(ageofdepth) > 0){
    depthseq <- sort(unique(c(ageofdepth, depthseq)))
  }

  # decide which models and point estimates should be used
  if(any(type==c(1, "int", "inter", "interp"))) type <- 1 else
  if(any(type==c(2, "reg", "regr", "poly", "polyn"))) type <- 2 else
  if(any(type==c(3, "spline", "spl"))) type <- 3 else
  if(any(type==c(4, "smooth", "sm"))) type <- 4 else
  if(any(type==c(5, "loess", "lowess"))) type <- 5 else
  if(any(type==c(6, "stineman"))) type <- 6
  
  if(est==1 || est==2) Est <- dat$mid1 else # 1 dummy, calculated later
  if(est==3) Est <- dat$mid1 else
  if(est==4) Est <- dat$wmn else
  if(est==5) Est <- dat$med else
  if(est==6) Est <- dat$mode else
  if(est==7) Est <- dat$mid2

  # remove outliers from the age-depth modelling
  if(length(outliers) > 0) {
    depths <- dat$depth[-outliers]
    errors <- dat$error[-outliers]
    calibs <- dat$calib[-outliers]
    Est <- Est[-outliers]
  } 
  else {
    depths <- dat$depth
    errors <- dat$error
    calibs <- dat$calib
  }

  # age-depth modelling with curves through sampled age estimates
  # in sections if one or more hiatuses are present
  if(length(hiatus) > 0) {
    
    allrange <- c(0,0,0,0)
    hiatusseq <- sort(c(range(depthseq), hiatus))
    
    for(i in 2:length(hiatusseq)) {
      
      cat(paste("\n section ", i-1, ",", sep=""))
      section <- depthseq[min(which(depthseq >= hiatusseq[i-1])) : max(which(depthseq <= hiatusseq[i]))]
      
      if(i>2) section <- section[-1]
      
      sel <- min(which(depths >= min(section))):max(which(depths <= max(section)))
      
      if(mixed.effect){
        if(length(outliers) > 0){
          smp <- .mixed.effect(its, depths, dat$cal[-outliers], dat$cage[-outliers], errors, calibs, Est, theta, f.mu, f.sigma, yrsteps, calibt)
        }
        else {
          smp <- .mixed.effect(its, depths, dat$cal, dat$cage, errors, calibs, est, theta, f.mu, f.sigma, yrsteps, calibt)
        }
      }
      else {
        smp <- .smpl(its, depths[sel], calibs[sel], Est[sel])
      }
      
      calrange <- .model.clam(type, smooth, its, wghts, depths[sel], errors[sel], section, prob, est, dat, smp, greyscale, remove.reverse, storedat, ageofdepth, BCAD)
      allrange <- rbind(allrange, calrange)
    }

    calrange <- allrange[2:nrow(allrange),]
  } 
  else {
    if(mixed.effect)
       if(length(outliers) > 0)
         smp <- .mixed.effect(its, depths, dat$cal[-outliers], dat$cage[-outliers], errors, calibs, est, theta, f.mu, f.sigma, yrsteps, calibt) else
           smp <- .mixed.effect(its, depths, dat$cal, dat$cage, errors, calibs, Est, theta, f.mu, f.sigma, yrsteps, calibt) else
             smp <- .smpl(its, depths, calibs, Est)
    calrange <- .model.clam(type, smooth, its, wghts, depths, errors, depthseq, prob, est, dat, smp, greyscale, remove.reverse, storedat, ageofdepth, BCAD)
  }
  dat$model <- approx(calrange[,1], (calrange[,2]+calrange[,3])/2, dat$depth)$y

  if(est==2){
    # Calibrated ages reduced to point estimates using midpoint of the calendar age.
    calrange[,4] <- (calrange[,2] + calrange[,3]) / 2
  }

  if(!BCAD && any(diff(calrange[,4]) < 0) || BCAD && any(diff(calrange[,4]) > 0)) {
    reversal <- TRUE
  }
  else{ 
    reversal <- FALSE
  }

  gfit <- round(.gfit(theta, f.mu, f.sigma, dat, calrange, outliers), 2)

  # re-correct the depths if slumps were applied
  if(length(slump) > 0)
  {
    dat <- .read.clam(name, ext, hpdsteps, yrsteps, prob, times, sep, BCAD, storedat, ignore, thickness, youngest, slump=c(), threshold, theta, f.mu, f.sigma, calibt, extradates, calcurve, bp) # read in the original dates again
    calrange <- calrange[which(calrange[,1] <= dmax),]
    d <- calrange[,1]
    for(i in 1:nrow(slump))
      {
        d[d > min(slump[i,])] <- d[d > min(slump[i,])] + (max(slump[i,]) - min(slump[i,]))
        dmax <- dmax + (max(slump[i,]) - min(slump[i,]))
        calrange[,1] <- d
        hiatus[hiatus > min(slump[i,])] <- hiatus[hiatus > min(slump[i,])] + (max(slump[i,]) - min(slump[i,]))
      }
  }

  # produce the age-depth plot, and a pdf copy if desired
  if(length(yrmin)==0) yrmin <- min(dat$mid1, calrange[,2])
  if(length(yrmax)==0) yrmax <- max(dat$mid1, calrange[,3])
  if(length(ageofdepth > 0)) layout(matrix(c(1,2,1,3), nrow=2), heights=c(.7,.3))
  .ageplot(yrmin, yrmax, dmin, dmax, revaxes, revd, revyr, yrlab, dlab, hiatus, depthseq, outliers, plotrange, BCAD, greyscale, if(length(greyscale)>0) chron else c(), C14col, outcol, outlsize, bestcol, rangecol, dat, calrange, depth, calhght, maxhght, mirror, calcol, slump, slumpcol, plotname, name, bty, mar, mgp, ash)

  # write files providing calibrated dates, age-model and settings
  colnames(calrange) <- c("Depth", paste("min.", 100*prob, "%range", sep=""), paste("max.", 100*prob, "%range", sep=""), "point")
  .write.clam(dat, runname, calrange, name, prob, type, remove.reverse, smooth, wghts, its, outliers, ignore, est, BCAD, yrsteps, every, decimals, accrate, depth, depthseq, hiatus, gfit, reversal, plotpdf, plotpng, yrmin, yrmax, dmin, dmax, dlab, yrlab, plotrange, greyscale, if(length(greyscale)>0) chron else c(), C14col, outcol, outlsize, bestcol, rangecol, calhght, maxhght, mirror, calcol, slump, slumpcol, revaxes, revyr, revd, calibt, youngest, extradates, plotname, calcurve, ccname, postbomb, pbnames, depths.file, bty, mar, mgp, ash)
  closeAllConnections()

  if(storedat)
  {
    calrange <<- calrange
    dat <<- dat
    smp <<- smp
  }

  # plot the age distribution of a provided depth
  if(length(ageofdepth) > 0)
  {
    if(revaxes)
      abline(v=ageofdepth, lty=2) else
        abline(h=ageofdepth, lty=2)
    xlim <- range(.ageofdepth)
    if(!BCAD) xlim <- xlim[2:1]
    hst <- density(.ageofdepth, n=max(1, max(xlim)-min(xlim)))
    yr <- seq(min(xlim), max(xlim), by=yrsteps)
    hst <- cbind(c(yr, max(xlim), min(xlim)), c(approx(hst$x, hst$y, yr)$y, 0, 0))
    plot(hst, type="n", main="", xlim=xlim, xlab=yrlab, ylab="")
    polygon(hst, col="grey")
    legend("topleft", paste(ageofdepth, depth), bty="n")
    layout(matrix(1))
    rng <- round(calrange[max(which(calrange[,1] <= ageofdepth)),])
    cat("\n  Age range of ", ageofdepth, " ", depth, ": ", rng[3], " to ", rng[2], ifelse(BCAD, " cal BC/AD", " cal BP"), " (", rng[3]-rng[2], " yr, ", prob, " % range)",  sep="")
  }

  # report the confidence ranges, the goodness-of-fit, and whether any age-reversals occurred
  rng <- round(calrange[,3]-calrange[,2])
  
  cat("\n  ", name, "'s ", 100*prob, "% confidence ranges span from ", min(rng), " to ", max(rng), " yr (average ", round(mean(rng)), " yr)", sep="")
  cat("\n  Fit (-log, lower is better):", gfit, "\n")
  
  if(reversal) cat("  Age reversals occurred. Try other model?\n")

  return(calrange)
}