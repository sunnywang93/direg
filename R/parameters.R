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
#' @export

theta_sheets <- function(X_list, tout, delta, e) {

  tout_minus <- tout |>
    apply(1, function(x) pmax(x - (delta/2 * e), 0)) |>
    t()

  tout_plus <- tout |>
    apply(1, function(x) pmin(x + (delta/2 * e), 1)) |>
    t()

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

  H_hat$e1[!is.finite(H_hat$e1)] <- 1

  list(
    H = pmin(pmax(H_hat, 0.1), 1),
    theta_e2 = theta_delta$e2,
    theta_e1 = theta_delta$e1
  )

}


#' Compute the Hölder regularity along the direction of input bases
#'
#' @param X_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param tout Matrix / dataframe with 2 dimensions, specifying all the
#' coordinates in the `x` and `y` axis.
#' @param delta Numeric, determining the spacings.
#' @param base_list List, with each element containing the coordinates of
#' each vector.
#' @returns List, containing a matrix with the regularity along each direction
#' of the respective bases functions.


H_sheets_dir <- function(X_list, tout, delta, base_list) {

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
                       ~(log(.x) - log(.y)) / (2 * log(2))
  ) |>
    purrr::map(~apply(.x, 2, function(x) pmin(pmax(x, 0.1), 1)))

  H_hat

}


