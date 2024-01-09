

#' Computation of regularity in the presence of noise
#'
#' Auxiliary function used in the computation of the regularity along some
#' directional vector.
#'
#' @param theta_2delta Numeric.
#' @param theta_delta Numeric.
#' @param sigma Numeric.
#' @returns Numeric.
#' @export

H_replace <- function(theta_2delta, theta_delta, sigma) {

  if(theta_2delta > 2*sigma**2 & theta_delta > 2*sigma**2) {
    (log(theta_2delta - 2*sigma**2) / (2*log(2))) -
      (log(theta_delta - 2*sigma**2) / (2*log(2)))
  } else {
    1
  }

}


identicalValue <- function(x,y) {
  if (identical(x,y)) {
      x
    } else {
      FALSE
    }
}





