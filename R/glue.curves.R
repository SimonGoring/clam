#' Bind calibration curves.
#' 
#' This function 'glues' two curves together, e.g., for calibrating southern hemisphere dates older than SHCal04
#' @param nh="IntCal09.14C"
#' @param sh="SHCal04.14C"
#' @param offset=c(56, 24)
#' @param name="gluedHemispheres.14C"
#' 
#' @return Outputs a table with a bound calibration curve.
#' @export
#' 
glue.curves <- function(nh="IntCal09.14C", sh="SHCal04.14C", offset=c(56, 24), name="gluedHemispheres.14C")
{
  nh <- nh[nh[,1] > max(sh[,1]),] # only use years beyond SHCal04
  nh[,2] <- nh[,2] + offset[1] # nh to sh, means
  nh[,3] <- sqrt(nh[,3]^2 + offset[2]^2) # errors, squared summed
  write.table(rbind(sh, nh), "gluedHemispheres.14C", sep="\t", row.names=F, col.names=F)
}
