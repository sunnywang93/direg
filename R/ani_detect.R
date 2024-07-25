#' Performs anisotropic detection by thresholding
#'
#' @param X_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **$X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param beta_n Numeric, number of random angles to take.
#' @param g_adj Numeric, the estimated and corrected function of angle.
#' @param alpha_adj Numeric, the estimated and correction angle.
#' @param H_min Numeric, the estimated minimum regularity.
#' @param H_max Numeric, the estimated maximum regularity.
#' @param delta Numeric, the spacings used for estimation.
#' @param sigma Numeric, the estimated noise.
#' @returns Numeric, where `0` indicates isotropy, and `1` indicates anisotropy.
#' @export

ani_detect <- function(X_list, beta_n, g_adj, alpha_adj,
                       H_min, H_max, delta, sigma,
                       H_pow = 1.1) {

  beta_vec <- runif(n = beta_n, min = alpha_adj + pi/4, max = alpha_adj + 3*pi/4)

  M0 <- length(X_list[[1]]$t)

  tsub <- seq(0, 1, length.out = floor(M0 / 3))

  tout <- expand.grid(t1 = tsub,
                      t2 = tsub)

  Hmin_beta <- purrr::map_dbl(beta_vec,
                              ~H_sheets_dir(X_list = X_list,
                                            tout = tout,
                                            delta = delta,
                                            base_list = list(c(cos(.x), sin(.x))),
                                            sigma = sigma))

  epsilon_beta <- mean(abs(outer(Hmin_beta, Hmin_beta, FUN = "-")))

  if(abs(H_max - H_min) > epsilon_beta + log(M0)^(-H_pow)) {
    message("Data is anisotropic")
    1
  } else {
    message("Data is isotropic")
    0
  }


}


