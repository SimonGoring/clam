#' Radiocarbon calibration on a student's t
#' 
#' See Christen and Perez 2009, Radiocarbon 51:1047-1059. Instead of assuming the standard Gaussian model (default in clam), a student t distribution can be used with two parameters. Christen and Perez 2009 suggest t.a = 3 and t.b = 4; this can be put as clam( calibt=c(3,4) )
#' 
#' @param t.a numeric.
#' @param t.b numeric.
#' @param f.cage numeric.
#' @param f.error numeric.
#' @param theta numeric.
#' @param f.mu numeric.
#' @param f.sigma numeric.
#' 
#' @return Probability distribution function.
#' @export
#'
.calibt <- function(t.a, t.b, f.cage, f.error, theta, f.mu, f.sigma){
  (t.b + ((f.cage-f.mu)^2) / (2*(f.sigma^2 + f.error^2))) ^ (-1*(t.a+0.5))
}