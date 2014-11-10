#'  Use student-t approach to model 14C calibration of a date and return a plot of the calibrated age.
#' 
#'  Using the caibration curve defined by the user, return the radiocarbon calibration curve (in cal yr BP) for a radiocarbon date.
#'  Using student's t to model the curve provides greater support to low probability dates at the margins
#'  of the calibration curve.
#' 
#' @param y numeric.  The radiocarbon year.
#' @param error numeric. The estimated error for the 14C age.
#' @param t.a numeric. The a parameter for the student's t distribution.
#' @param t.b numeric. The b parameter for the student's t distribution.
#' @param cc numeric. Should this value be calibrated, or is it already calibrated?
#' @param postbomb numeric.  A value of 1 - 5 related to the postbomb curve that should be used.  See details.
#' @param cc1 character. One of IntCal13, Marine13, SHCal13, or ConstCal.
#' @param cc2 character. As above.
#' @param cc3 character. As above.
#' @param cc4 character. As above.
#' @param Cutoff numeric.
#' 
#' @return A plot of the calibrated age.
#' @export

student.t <- function(y=2450, error=50, t.a=3, t.b=4, cc=1, postbomb=c(), cc1="IntCal13",
                      cc2="Marine13", cc3="SHCal13", cc4="ConstCal", Cutoff=1e-5)
{
  if(cc==0)
  {
    cc <- seq(y-10*error, y+10*error, length=1e3)
    cc <- cbind(cc, cc, rep(0, length(cc)))
  } else
  {
    if(cc1=="IntCal13") cc1 <- read.table("data/IntCal13.14C") else
      cc1 <- read.csv(cc1)[,1:3]
    if(cc2=="Marine13") cc2 <- read.table("data/Marine13.14C") else
      cc2 <- read.csv(cc2)[,1:3]
    if(cc3=="SHCal13") cc3 <- read.table("data/SHCal13.14C") else
      cc3 <- read.table(cc3)[,1:3]
    if(cc4 != "ConstCal") cc4 <- read.table(paste('data/', cc4, sep = ''))[,1:3]
    if(cc==1) cc <- cc1 else if(cc==2) cc <- cc2 else if(cc==3) cc <- cc3 else cc <- cc4
  }
  
  if(y < 0)
    if(length(postbomb) > 0)
    {
      if(postbomb==1) bomb <- read.table("postbomb_NH1.14C")[,1:3] else
        if(postbomb==2) bomb <- read.table("postbomb_NH2.14C")[,1:3] else
          if(postbomb==3) bomb <- read.table("postbomb_NH3.14C")[,1:3] else
            if(postbomb==4) bomb <- read.table("postbomb_SH1-2.14C")[,1:3] else
              if(postbomb==5) bomb <- read.table("postbomb_SH3.14C")[,1:3] else
                stop("Warning, cannot find postbomb curve #", postbomb, " (use values of 1 to 5 only)")
      bomb.x <- seq(max(bomb[,1]), min(bomb[,1]), by=-.1) # interpolate
      bomb.y <- approx(bomb[,1], bomb[,2], bomb.x)$y
      bomb.z <- approx(bomb[,1], bomb[,3], bomb.x)$y
      bomb <- cbind(bomb.x, bomb.y, bomb.z, deparse.level=0)
      if(info$postbomb < 4)
        cc1 <- rbind(bomb, cc1, deparse.level=0) else
          cc3 <- rbind(bomb, cc3, deparse.level=0)
    }
  
  norm.cal <- dnorm(cc[,2], y, sqrt(cc[,3]^2+error^2))
  t.cal <- (t.b + ((y-cc[,2])^2) / (2*(cc[,3]^2 + error^2))) ^ (-1*(t.a+0.5))
  
  norm.cal <- cbind(cc[,1], norm.cal/sum(norm.cal))
  acc <- which(norm.cal[,2] >= Cutoff)
  if(length(acc) > 1) norm.cal <- norm.cal[acc,]
  
  t.cal <- cbind(cc[,1], t.cal/sum(t.cal))
  acc <- which(t.cal[,2] >= Cutoff)
  if(length(acc) > 1) t.cal <- t.cal[acc,]
  t.cal <- cbind(c(min(t.cal[,1]), t.cal[,1], max(t.cal[,1])), c(0, t.cal[,2], 0))
  
  plot(norm.cal, type="l", xlab="cal BP", xlim=range(c(t.cal[,1], norm.cal[,1]))[2:1], ylab="", ylim=c(0, max(t.cal[,2], norm.cal[,2])), col=2, lwd=1.5)
  polygon(t.cal, col=rgb(0,0,0,.25), border=rgb(0,0,0,.5))
  legend("topright", "Gaussian", text.col=2, bty="n")
  legend("topright", paste("\nstudent-t (a=", t.a, ", b=", t.b, ")", sep=""), bty="n", text.col=grey(.4))
}
