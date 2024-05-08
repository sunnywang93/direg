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
#' the regularity.
#' @param xout Vector, containing the evaluation points along one 1 dimension
#' to compute the regularity. The cartesian product is taken to produce the
#' 2D grid.
#' @returns List, containing the identified angle and the associated estimated
#' regularity averaged over the grid of spacings.
#' @export

identify_angle <- function(angles, X_list, dout, xout) {

  # Construct the two basis vectors
  v1_cot <- c(cos(angles$alpha_acot), sin(angles$alpha_acot))
  v1_tan <- c(cos(angles$alpha_atan), sin(angles$alpha_atan))
  # Construct the two reflected basis vectors
  v1_cot_ref <- c(cos(pi - angles$alpha_acot),
                  sin(pi - angles$alpha_acot))
  v1_tan_ref <- c(cos(pi - angles$alpha_atan),
                  sin(pi - angles$alpha_atan))
  # Compute the regularity along each basis vector along a grid of deltas
  noise <- estimate_sigma(X_list = X_list)

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

  # Compute the average of regularities
  H_avg <- max(Reduce('+', H_v) / length(H_v))

  # Return the unique angle that maximises the regularity
  if(mode_idx == 1) {
    alpha <- angles$alpha_acot
    names(alpha) <- "alpha_acot"
  } else if(mode_idx == 2) {
    alpha <- pi - angles$alpha_acot
    names(alpha) <- "alpha_acot_ref"
  } else if(mode_idx == 3) {
    alpha <- angles$alpha_atan
    names(alpha) <- "alpha_atan"
  } else {
    alpha <- pi - angles$alpha_atan
    names(alpha) <- "alpha_atan_ref"
  }


  list(alpha = alpha,
       H_max = H_avg)

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
                    na.rm = TRUE)**(1 / (2 * H_list$H))

  list(g_hat = tan_alpha,
      alpha_acot = pracma::acot(tan_alpha),
      alpha_atan = atan(tan_alpha),
      H_min = H_list$H)

}


#' Compute the correction for estimating equation of angle
#'
#' Correction is based on the remainder term that arises from the estimating
#' equation used for the directional regularity. Should only be used once to
#' avoid introducing additional dependence between estimates.
#'
#' @param g_hat Numeric, the preliminary estimate of the tangent or cotangent of the angle.
#' @param alpha_hat Numeric, the preliminary estimate of the identified angle.
#' @param delta Numeric, the spacing used for estimation.
#' @param Hmax Numeric, the maximum regularity.
#' @param Hmin Numeric, the minimum regularity.
#' @returns Numeric, the corrected angle.
#' @export

angle_correct <- function(g_hat, alpha_hat, delta, Hmax, Hmin) {

  if(names(alpha_hat) == "alpha_acot" | names(alpha_hat) == "alpha_acot_ref") {
    f_alpha_denom <- abs(sin(alpha_hat) * delta)^(2*H_min) +
      abs(cos(alpha_hat) * delta)^(2*H_max)
    f_alpha_1 <- abs(sin(alpha_hat) * delta)^(2*H_min) / f_alpha_denom
    f_alpha_2 <- (abs(tan(alpha_hat))^(2*H_min) *
                    abs(sin(alpha_hat*delta))^(2*H_max)) / f_alpha_denom
    f_alpha <- (f_alpha_1 + f_alpha_2)^(1/(2*H_min))
  } else {
    f_alpha_denom <- abs(sin(alpha_hat) * delta)^(2*H_max) +
      abs(cos(alpha_hat) * delta)^(2*H_min)
    f_alpha_1 <- abs(cos(alpha_hat) * delta)^(2*H_min) / f_alpha_denom
    f_alpha_2 <- (abs(pracma::cot(alpha_hat))^(2*H_min) *
                    abs(cos(alpha_hat) * delta)^(2*H_max)) / f_alpha_denom
    f_alpha <- (f_alpha_1 + f_alpha_2)^(2*H_min)
  }

  alpha_adj <- pracma::acot(g_hat / f_alpha)

  # Use previous identification before correction to obtain the right alpha
  if(names(alpha_hat) == "alpha_acot_ref" | names(alpha_hat) == "alpha_atan_ref") {
    alpha_adj <- pi - alpha_adj
  }

  alpha_adj

}



#' Performs a change-of-basis for bivariate functional data
#'
#' Given a rotational matrix, a change of basis is performed for the sampling
#' points of functional data in the form of bivariate sheets. Relevant for
#' anisotropic functional data with common design.
#'
#' @param Y_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **$X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param rot_mat Matrix, containing the coordinates of the linear transformation.
#' Typically a rotation matrix.
#' @returns List, containing the transformed sampling points and original
#' observations.
#' @export

change_basis <- function(Y_list, rot_mat) {

  tout <- as.matrix(
    expand.grid(t1 = Y_list[[1]]$t, t2 = Y_list[[1]]$t)
  )

  Ralpha_t <- t(apply(tout, 1, function(tm) crossprod(t(rot_mat), tm)))


  purrr::map(Y_list,
             ~list(t = Ralpha_t,
                   X = .x$X))

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


