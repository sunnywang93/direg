
#' Estimates the angle between two basis vectors
#'
#' @param X_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **$X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param xout Vector, containing the evaluation points at which
#' the angle should be computed.
#' @param delta Numeric, determining the spacings.
#' @export

estimate_angle <- function(X_list, xout, delta) {

  tout <- expand.grid(xout, xout)

  H_list <- H_sheets(X_list = X_list,
                     tout = tout,
                     delta = delta)


  tan_alpha <- (H_list$theta_e2 / H_list$theta_e1)**(1 / (2 * mean(H_list$H)))

  list(alpha_cot = pracma::acot(tan_alpha),
       alpha_tan = atan(tan_alpha))

}



