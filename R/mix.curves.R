#' Create a new calibration curve if dates are of mixed source.
#' 
#' If two curves need to be 'mixed' to calibrate, then re-generate a calibration curve. e.g. for dates of mixed terrestrial and marine carbon sources
#' 
#' @param ratio=.5
#' @param cc1="IntCal13.14C"
#' @param cc2="Marine13.14C"
#' @param name="mixed.14C"
#' @param offset=c(0,0)
#' 
#' @return Output to a new table at location defined in \code{name}.
#' @export
#' 
mix.curves <- function(ratio=.5, cc1="IntCal13.14C", cc2="Marine13.14C", name="mixed.14C", offset=c(0,0))
{
  cc1 <- read.table(cc1)
  cc2 <- read.table(cc2)
  cc2.mu <- approx(cc2[,1], cc2[,2], cc1[,1], rule=2)$y + offset[1] # interpolate cc2 to the calendar years of cc1
  cc2.error <- approx(cc2[,1], cc2[,3], cc1[,1], rule=2)$y
  cc2.error <- sqrt(cc2.error^2 + offset[2]^2)
  mu <- ratio * cc1[,2] + (1-ratio) * cc2.mu
  error <- ratio * cc1[,3] + (1-ratio) * cc2.error
  write.table(cbind(cc1[,1], mu, error), name, row.names=FALSE, col.names=FALSE, sep="\t")
}
