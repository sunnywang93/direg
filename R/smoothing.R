
#' Calculate plug-in bandwidth for adaptive bivariate kernel smoothing
#'
#' Currently implemented for surfaces with common design sampling points.
#'
#' @param H1 Numeric, containing the Hölder regularity along the first dimension.
#' @param H2 Numeric, containing the Hölder regularity along the second dimension.
#' @param L1 Numeric, containing the Hölder constant along the first dimension.
#' @param L2 Numeric, containing the Hölder constant along the second dimension.
#' @param sigma Numeric, containing the noise level.
#' @param k Numeric, constant for the lower bound of the kernel.
#' @param M0 Numeric, number of points along each curve.
#' @returns Vector, containing the estimated bandwidths along each dimension.
#' @export

bw_smooth <- function(H1, H2, L1, L2, sigma, k, M0) {

  Lambda1_k1 <- ((k^2 * sigma^2) / (4 * L1 * H1))**(2*H2 + 1)

  Lambda1_k2 <- (4 * L2 * H2) / (k^2 * sigma^2)

  Lambda1 <- (Lambda1_k1 * Lambda1_k2)**(1 / (4*H1*H2 + 2*H1 + 2*H2))

  rate1 <- M0**(-H2 / (2*H2*H1 + H1 + H2))

  h1 <- Lambda1 * rate1

  Lambda2_k1 <- ((k^2 * sigma^2) / (4 * L2 * H2))**(2*H1 + 1)

  Lambda2_k2 <- (4 * L1 * H1) / (k^2 * sigma^2)

  Lambda2 <- (Lambda2_k1 * Lambda2_k2)**(1 / (4*H1*H2 + 2*H1 + 2*H2))

  rate2 <- M0**(-H1 / (2*H2*H1 + H1 + H2))

  h2 <- Lambda2 * rate2

  c(h1 = h1, h2 = h2)


}

#' Bi-variate kernel smoothing for functional data on common design
#'
#' Multiplicative kernels are used, with a diagonal bandwidth matrix.
#' Warning! If smoothing is desired after a change of basis from considering
#' the directional regularity, remember to also change the evaluation points
#' in addition to the input list of sampling and observed points.
#'
#' @param Y_list List, containing the following elements:
#' - **$t** Vector of sampling points.
#' - **$X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param bw_vec Vector of bandwidths, corresponding to the diagonal elements of
#' the bandwidth matrix.
#' @param xout Vector of evaluation points along one dimension. Smoothing is performed
#' on the cartesian product of `xout` with itself.
#' @returns List, corresponding to the smooth curves.
#' @export


ksmooth_bi <- function(Y_list, bw_vec, xout) {

  tobs <- expand.grid(t1 = Y_list[[1]]$t, t2 = Y_list[[1]]$t)

  tout <- expand.grid(t1 = xout, t2 = xout)

  # Faster than outer
  weights <- vapply(seq_along(tobs$t1),
                    function(x) epa_kernel((tobs$t1[x] - tout$t1) / bw_vec[1])
                    * epa_kernel((tobs$t2[x] - tout$t2) / bw_vec[2]),
                    FUN.VALUE = numeric(length(tout$t1))
                    ) |>
    (\(x) t(sweep(x, MARGIN = 1, STATS = rowSums(x), FUN = "/")))() |>
    (\(x) replace(x, is.nan(x), 0))()

  purrr::map(Y_list,
             ~matrix(crossprod(c(.x$X), weights),
                     nrow = length(xout),
                     ncol = length(xout)
                     )
             )

}



