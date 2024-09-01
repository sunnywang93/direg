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
#' @param delta_grid Vector, the grid of spacings used for estimation.
#' @param sigma Numeric, the estimated noise.
#' @param H_pow Numeric, the power of the rate term.
#' @returns List, containing the following elements:
#'  - **$indic**: Numeric, with 1 indicator anisotropy and 0 otherwise.
#'  - **$tau** : Numeric, the threshold value used to test for anisotropy.
#' @export

ani_detect <- function(X_list, beta_n, g_adj, alpha_adj,
                       H_min, H_max, delta_grid, sigma,
                       H_pow = 1/3) {

  beta_vec <- runif(n = beta_n,
                    min = alpha_adj + pi/4,
                    max = alpha_adj + 3*pi/4)

  M0 <- length(X_list[[1]]$t)^2

  tsub <- seq(0, 1, length.out = floor(sqrt(M0) / 3))

  tout <- expand.grid(t1 = tsub,
                      t2 = tsub)

  base_beta <- purrr::map2(cos(beta_vec), sin(beta_vec),
                           ~c(.x, .y))
  base_beta_ortho <- purrr::map2(cos(beta_vec + pi/2), sin(beta_vec + pi/2),
                                 ~c(.x, .y))

  Hmin_beta <- purrr::map(delta_grid,
                          ~H_sheets_dir(X_list = X_list,
                                        tout = tout,
                                        delta = .x,
                                        base_list = base_beta,
                                        sigma = sigma)) |>
    (\(x) Reduce('+', x) / length(x))()

  Hmin_beta_ortho <- purrr::map(delta_grid,
                                ~H_sheets_dir(X_list = X_list,
                                              tout = tout,
                                              delta = .x,
                                              base_list = base_beta_ortho,
                                              sigma = sigma)) |>
    (\(x) Reduce('+', x) / length(x))()

  epsilon_beta <- mean(abs(Hmin_beta - Hmin_beta_ortho), na.rm = TRUE)
  tau <- epsilon_beta + exp(-log(M0)^H_pow)

  if(abs(H_max - H_min) > tau) {
    #message("Data is anisotropic")
    indic <- 1
  } else {
    #message("Data is isotropic")
    indic <- 0
  }

  list('indic' = indic,
       'tau' = tau)


}






