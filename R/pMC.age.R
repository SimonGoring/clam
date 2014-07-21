#' Calculate 14C ages from pmC values
#' 
#' @param mn
#' @param sdev
#' @param ratio=100
#' @param decimals=0
#' 
#' @return Rounded estimates of the 14C dates.
#' @export
#' 
pMC.age <- function(mn, sdev, ratio=100, decimals=0)
{
  y <- -8033*log(mn/ratio)
  sdev <- y - -8033*log((mn+sdev)/ratio)
  round(c(y, sdev), decimals)
}

