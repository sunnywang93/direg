#' Computes the mean squared deviations of a surface on a grid on points
#'
#'
#' @param X_list List, containing the following elements:
#' -**$t** Vector of sampling points,
#' -**X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param tout Data frame of points on which the thetas should be computed,
#' resulting from `expand.grid`.
#' @param delta Numeric, determining the spacings.
#' @param e Vector, containing the coordinates of the directional vector.
#' @returns Matrix, containing the mean squared deviations on the evaluation
#' points `tout`.
#' @export

theta_sheets <- function(X_list, tout, delta, e) {

  tout_minus <- cbind(pmin(pmax(tout[, 1] - (delta/2 * e[1]), 0), 1),
                      pmin(pmax(tout[, 2] - (delta/2 * e[2]), 0), 1)
                      )

  tout_plus <- cbind(pmax(pmin(tout[, 1] + (delta/2 * e[1]), 1), 0),
                     pmax(pmin(tout[, 2] + (delta/2 * e[2]), 1), 0)
                     )

  X_minus <- purrr::map(X_list,
                        ~pracma::interp2(x = .x$t,
                                         y = .x$t,
                                         Z = .x$X,
                                         xp = tout_minus[, 2],
                                         yp = tout_minus[, 1],
                                         method = "nearest") |>
                          matrix(nrow = sqrt(nrow(tout)),
                                 ncol = sqrt(nrow(tout))))

  X_plus <- purrr::map(X_list,
                       ~pracma::interp2(x = .x$t,
                                        y = .x$t,
                                        Z = .x$X,
                                        xp = tout_plus[, 2],
                                        yp = tout_plus[, 1],
                                        method = "nearest") |>
                         matrix(nrow = sqrt(nrow(tout)),
                                ncol = sqrt(nrow(tout))))

  purrr::map2(X_minus, X_plus,
              ~(.x - .y)**2) |>
    (\(x) Reduce('+', x) / length(x))()

}

#' Computes the Hölder regularity of bivariate functional data
#'
#' @param X_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param tout Matrix / dataframe with 2 dimensions, specifying all the
#' coordinates in the `x` and `y` axis.
#' @param delta Numeric, determining the spacings.
#' @returns List, containing the elements:
#' - **H** Matrix, containing the estimated regularity
#' on the dimensional grid `tout`.
#' - **theta_e2** Matrix, containing the thetas on the `e2` basis computed with
#' delta.
#' - **theta_e1** Matrix, containing the thetas on the `e1` basis computed with
#' delta.
#' @export

H_sheets <- function(X_list, tout, delta) {

  e <- list(e1 = c(1, 0),
            e2 = c(0, 1))

  theta_2delta <- purrr::map(e,
                             ~theta_sheets(X_list = X_list,
                                           tout = tout,
                                           delta = delta * 2,
                                           .x))

  theta_delta <- purrr::map(e,
                            ~theta_sheets(X_list = X_list,
                                          tout = tout,
                                          delta = delta,
                                          .x))

  H_hat <- purrr::map2(theta_2delta, theta_delta,
                       ~(log(.x) - log(.y)) / (2 * log(2))
                       ) |>
    (\(x) Reduce(pmin, x))()

  list(
    H = mean(pmin(pmax(H_hat, 0.1), 1)),
    theta_e2 = theta_delta$e2,
    theta_e1 = theta_delta$e1
  )

}


#' Compute the Hölder regularity along the direction of input bases
#'
#' @param X_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **$X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param tout Matrix / dataframe with 2 dimensions, specifying all the
#' coordinates in the `x` and `y` axis.
#' @param delta Numeric, determining the spacings.
#' @param base_list List, with each element containing the coordinates of
#' each vector.
#' @param sigma Numeric, indicating the noise level, for example from the
#' `estimated_sigma` function.
#' @returns Vector, containing the estimated regularity along each direction
#' of `base_list`. The order of outputs corresponds to the order of inputs in
#' `base_list`.
#' @export


H_sheets_dir <- function(X_list, tout, delta, base_list, sigma) {


  theta_2delta <- purrr::map(base_list,
                             ~theta_sheets(X_list = X_list,
                                           tout = tout,
                                           delta = delta * 2,
                                           .x))

  theta_delta <- purrr::map(base_list,
                            ~theta_sheets(X_list = X_list,
                                          tout = tout,
                                          delta = delta,
                                          .x))

  H_hat <- purrr::map2(theta_2delta, theta_delta,
                       ~matrix(mapply(H_replace, .x, .y, sigma),
                               nrow = nrow(.x),
                               ncol = ncol(.x)
                               )
                       )

  purrr::map_dbl(H_hat, ~mean(pmax(.x, 0.1), na.rm = TRUE))

}


#' Estimate the noise at the observed points of functional surfaces
#'
#' @param X_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @returns Matrix, where the number of rows and columns are each given by the
#' length of the sampling points.
#' @export

estimate_sigma <- function(X_list) {

  # Identify nearest neighbours for each observed point
  tout <- expand.grid(t1 = X_list[[1]]$t, t2 = X_list[[1]]$t)
  min_idx <- sapply(seq_len(nrow(tout)), function(y) {
    idx <- which.min(
      abs(tout$t1[y] - tout$t1[-y]) + abs(tout$t2[y] - tout$t2[-y])
      )
    ifelse(idx >= y, idx+1, idx)
    })

  # Compute noise by averaging over all curves
  lapply(X_list, function(curve) {
    sapply(seq_along(curve$X), function(id) {
      (curve$X[id] - curve$X[min_idx[id]])**2
    })
  }) |>
    (\(x) mean(sqrt(Reduce('+', x) / (2*length(x)))))()

}


#' Estimates the Hölder constant along some directional vector
#'
#' Computation is done using the mean-squared continuity that arises from the class
#' of stochastic process considered. Estimates are first computed on a grid of time
#' points, before the average over the time points are taken to stabilise estimates.
#'
#' @param X_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **$X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param xout Vector, containing the evaluation points at which
#' the estimates should be computed along one dimension. Computation is done
#' over the cartesian product of `xout` with itself.
#' @param delta Numeric, determining the spacings.
#' @param v Vector, determining the direction at which the Hölder constant should
#' be estimated.
#' @param Hv Numeric, containing the Hölder regularity along direction `v`.
#' @returns Numeric, containing the Hölder constant.
#' @export

L_sheets <- function(X_list, xout, delta, v, Hv) {

  tout <- expand.grid(t1 = xout, t2 = xout)

  theta_delta <- theta_sheets(X_list = X_list,
                              tout = tout,
                              delta = delta,
                              v)

  mean(theta_delta / delta**(2*Hv))

}




