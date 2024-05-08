

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

#' Computes the Epanechnikov Kernel
#'
#' @param y Vector of points.
#' @returns Vector.
#' @export

epa_kernel <- function(y) {
  3/4 * (1 - y^2) * (abs(y) <= 1)
}


#' Computes the risk of the true and estimated tangent or cotangent of alpha,
#' and their relative risk.
#'
#' @param alpha Numeric, the true angle.
#' @param H1 Numeric, the regularity along the first basis.
#' @param H2 Numeric, the regularity along the second basis.
#' @param delta Numeric, the spacings used in estimation.
#' @param ghat Numeric, the estimated function of the angle.
#' @returns Numeric, the true tangent or cotangent of the angle.
#' @export

g_risk <- function(alpha, H1, H2, delta, ghat) {

  g_alpha_num <- abs(sin(alpha) * delta)^(2*H1) +
    abs(cos(alpha) * delta)^(2*H2)

  g_alpha_denom <- abs(cos(alpha))^(2*H1) * delta^(2*H1) +
    abs(sin(alpha))^(2*H2) * delta^(2*H2)

  g_true <- (g_alpha_num / g_alpha_denom)^(1/ (2*min(H1, H2)))

  if(H1 > H2) {
    ghat_risk <- abs(ghat - pracma::cot(alpha))
    g_risk <- abs(g_true - pracma::cot(alpha))
  } else {
    ghat_risk <- abs(ghat - tan(alpha))
    g_risk <- abs(g_true - tan(alpha))
  }

  list(
    g_true = g_true,
    ghat_risk = ghat_risk,
    g_risk = g_risk,
    g_rel = ghat_risk / g_risk
  )

}




