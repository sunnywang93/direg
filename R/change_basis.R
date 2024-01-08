#' Find the unique angle between the two basis vectors that provides the
#' maximising regularity
#'
#' @param angles Vector, containing the arccot and arctan of the angle, for
#' example outputted by the function `estimate_angle`.
#' @param X_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **$X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param dout Vector, containing the grid of spacings (i.e delta) to compute
#' the regularity. Plurality vote over the grid to find the maximising angles.
#' @param xout Vector, containing the evaluation points along one 1 dimension
#' to compute the regularity. The cartesian product is taken to produce the
#' 2D grid.
#' @returns Numeric, containing the identified angle.
#' @export

identify_angle <- function(angles, X_list, dout, xout) {

  # Construct the two basis vectors
  v1_cot <- c(cos(angles["alpha_cot"]), sin(angles["alpha_cot"]))
  v1_tan <- c(cos(angles["alpha_tan"]), sin(angles["alpha_tan"]))
  # Construct the two reflected basis vectors
  v1_cot_ref <- c(cos(pi - angles["alpha_cot"]),
                  sin(pi - angles["alpha_cot"]))
  v1_tan_ref <- c(cos(pi - angles["alpha_tan"]),
                  sin(pi - angles["alpha_tan"]))
  # Compute the regularity along each basis vector along a grid of deltas
  noise <- mean(estimate_sigma(X_list = X_list))

  tout <- expand.grid(t1 = xout, t2 = xout)

  H_v <- purrr::map(dout,
                    ~H_sheets_dir(X_list = X_list,
                     tout = tout,
                     delta = .x,
                     base_list = list(v1_cot, v1_cot_ref,
                                      v1_tan, v1_tan_ref),
                     sigma = noise))
  # Compute the sum of regularities across the grid
  mode_idx <- which.max(Reduce('+', H_v))

  # Return the unique angle that maximises the regularity
  if(mode_idx == 1) {
    angles["alpha_cot"]
  } else if(mode_idx == 2) {
    pi - angles["alpha_cot"]
  } else if(mode_idx == 3) {
    angles["alpha_tan"]
  } else {
    pi - angles["alpha_tan"]
  }

}



#' Estimates the angle between two basis vectors up to a reflection
#'
#' @param X_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **$X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param xout Vector, containing the evaluation points at which
#' the angle should be computed.
#' @param delta Numeric, determining the spacings.
#' @returns Vector, containing the angles identified by arccot and arctan.
#' @export

estimate_angle <- function(X_list, xout, delta) {

  tout <- expand.grid(xout, xout)

  H_list <- H_sheets(X_list = X_list,
                     tout = tout,
                     delta = delta)

  tan_alpha <- mean(H_list$theta_e2 / H_list$theta_e1,
                    na.rm = TRUE)**(1 / (2 * mean(H_list$H, na.rm = TRUE)))

  c(alpha_cot = pracma::acot(tan_alpha),
    alpha_tan = atan(tan_alpha))

}


alpha_avg <- function(X, npart) {

  xy <- expand.grid(x = seq_len(nrow(X)), y = seq_len(ncol(X)))

  breaks_xy <- sapply(seq_len(npart), function(i) ceiling(i * (nrow(X)/npart)))

  cuts_x <- cut(xy$x,
                breaks = c(1, breaks_xy),
                labels = FALSE,
                include.lowest = TRUE)
  cuts_y <- cut(xy$y,
                breaks = c(1, breaks_xy),
                labels = FALSE,
                include.lowest = TRUE)

  index_matrix <- matrix(cuts_x + (cuts_y - 1) * max(cuts_x),
                         nrow = nrow(X),
                         byrow = TRUE)

  ave(X, index_matrix, FUN = mean)

}


