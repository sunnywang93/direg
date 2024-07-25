
#' Calculate plug-in bandwidth for adaptive bivariate kernel smoothing
#'
#' Currently implemented for surfaces with common design sampling points.
#'
#' @param H1 Numeric, containing the Hölder regularity along the first dimension.
#' @param H2 Numeric, containing the Hölder regularity along the second dimension.
#' @param M0 Numeric, number of points along each curve.
#' @returns Vector, containing the estimated bandwidths along each dimension.
#' @export

bw_smooth <- function(H1, H2, M0, rate = TRUE) {

  rate1 <- M0**(-H2 / (2*H2*H1 + H1 + H2))
  rate2 <- M0**(-H1 / (2*H2*H1 + H1 + H2))
  h1 <- rate1
  h2 <- rate2
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
#' @param xout Matrix of evaluation points in two dimensions.
#' @returns List, corresponding to the smooth curves.
#' @export


ksmooth_bi <- function(Y_list, bw_vec, tobs, tout) {


  # Faster than outer
  weights <- vapply(seq_along(tobs[, 1]),
                    function(x) epa_kernel((tobs[, 1][x] - tout[, 1]) / bw_vec[1])
                    * epa_kernel((tobs[, 2][x] - tout[, 2]) / bw_vec[2]),
                    FUN.VALUE = numeric(length(tout[, 1]))
                    ) |>
    (\(x) t(sweep(x, MARGIN = 1, STATS = rowSums(x), FUN = "/")))() |>
    (\(x) replace(x, is.nan(x), 0))()

  purrr::map(Y_list,
             ~matrix(crossprod(c(.x$X), weights),
                     nrow = sqrt(length(tout[, 1])),
                     ncol = sqrt(length(tout[, 1]))
                     )
             )

}



# ksmooth_test <- function(Y_list, bw_vec, tobs, tout) {
#
#
#   # Faster than outer
#   weights <- vapply(seq_len(nrow(tobs)),
#                     function(x) multi_epa(
#                       cbind(((tobs[, 1][x] - tout[, 1])^2 / bw_vec[1]),
#                             ((tobs[, 2][x] - tout[, 2])^2 / bw_vec[2])
#                       )
#                     ),
#                     FUN.VALUE = numeric(length(tout[, 1]))
#   ) |>
#     (\(x) t(sweep(x, MARGIN = 1, STATS = rowSums(x), FUN = "/")))() |>
#     (\(x) replace(x, is.nan(x), 0))()
#
#   purrr::map(Y_list,
#              ~matrix(crossprod(c(.x$X), weights),
#                      nrow = sqrt(length(tout[, 1])),
#                      ncol = sqrt(length(tout[, 1]))
#              )
#   )
#
# }
#
#
# multi_epa <- function(Y) {
#   apply(Y, 1, function(x) 3/4 * (1 - x[1]^2 - x[2]^2) * (x[1]^2 + x[2]^2 <= 1))
# }
#
#
#
#
