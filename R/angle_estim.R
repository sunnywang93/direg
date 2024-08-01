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
#' @param sigma Numeric, the standard deviation of the noise,
#' for example outputted by `estimate_angle`.
#' @returns List, containing the identified angle and the associated estimated
#' regularity averaged over the grid of spacings.
#' @export

identify_angle <- function(angles, X_list, dout, xout, sigma) {

  # Construct the two basis vectors
  v1_cot <- c(cos(angles$alpha_acot), sin(angles$alpha_acot))
  v1_tan <- c(cos(angles$alpha_atan), sin(angles$alpha_atan))
  # Construct the two reflected basis vectors
  v1_cot_ref <- c(cos(pi - angles$alpha_acot),
                  sin(pi - angles$alpha_acot))
  v1_tan_ref <- c(cos(pi - angles$alpha_atan),
                  sin(pi - angles$alpha_atan))
  # Compute the regularity along each basis vector along a grid of deltas
  tout <- expand.grid(t1 = xout, t2 = xout)

  H_v <- purrr::map(dout,
                    ~H_sheets_dir(X_list = X_list,
                     tout = tout,
                     delta = .x,
                     base_list = list(v1_cot, v1_cot_ref,
                                      v1_tan, v1_tan_ref),
                     sigma = sigma))
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

  theta_ratio_num <- H_list$theta_e2
  theta_ratio_denom <- H_list$theta_e1

  theta_ratio <- theta_ratio_num / theta_ratio_denom

  g_alpha <- theta_ratio^(1 / (2 * H_list$H))


  list(g_hat = g_alpha,
      alpha_acot = pracma::acot(g_alpha),
      alpha_atan = atan(g_alpha),
      H_min = H_list$H,
      sigma_hat = H_list$sigma_hat)

}


#' Compute the correction for estimating equation of angle
#'
#' Correction is based on the remainder term that arises from the estimating
#' equation used for the directional regularity. Should only be used once to
#' avoid introducing additional dependence between estimates.
#'
#' @param X_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **$X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param xout Vector, containing the evaluation points at which the angle should be computed.
#' @param g_hat Numeric, the preliminary estimate of the tangent or cotangent of the angle.
#' @param alpha_hat Numeric, the preliminary estimate of the identified angle.
#' @param delta Numeric, the spacing used for estimation.
#' @param Hmax Numeric, the maximum regularity.
#' @param Hmin Numeric, the minimum regularity.
#' @param sigma Numeric, the standard deviation of the noise,
#' for example outputted by `estimate_angle`.
#' @returns List, containing the corrected angle, corrected function of angle and
#' the correction term.
#' @export

angle_correct <- function(X_list, xout, g_hat, alpha_hat, delta, Hmax, Hmin, sigma) {

  tout <- expand.grid(xout, xout)

  base_list <- list(u1 = c(cos(alpha_hat), sin(alpha_hat)),
                    u2 = c(-sin(alpha_hat), cos(alpha_hat)))

  theta_u <- purrr::map(base_list,
                        ~mean(theta_sheets(X_list = X_list,
                                           tout = tout,
                                           delta = delta,
                                           e = .x,
                                           sigma),
                              na.rm = TRUE)
                        )


  if(names(alpha_hat) == "alpha_acot" | names(alpha_hat) == "alpha_acot_ref") {

    f_alpha_num <- (abs(sin(alpha_hat))^(2*Hmax) / abs(cos(alpha_hat))^(2*Hmin)) *
      (theta_u$u1 / theta_u$u2)

    f_alpha_denom <- (abs(cos(alpha_hat))^(2*Hmax) / abs(sin(alpha_hat))^(2*Hmin)) *
      (theta_u$u1 / theta_u$u2)

    f_alpha <- ((1 + f_alpha_num) / (1 + f_alpha_denom))^(1 / (2 * Hmin))

  } else {

    f_alpha_num <- (abs(cos(alpha_hat))^(2*Hmax) / abs(sin(alpha_hat))^(2*Hmin)) *
      delta^(2*(Hmax - Hmin))

    f_alpha_denom <- (abs(sin(alpha_hat))^(2*Hmax) / abs(cos(alpha_hat))^(2*Hmin)) *
      delta^(2*(Hmax - Hmin))

    f_alpha <- ((1 + f_alpha_num) / (1 + f_alpha_denom))^(1 / 2 * Hmin)

  }

  g_adj <- g_hat / f_alpha

  if(names(alpha_hat) == "alpha_acot" | names(alpha_hat) == "alpha_acot_ref") {
    alpha_adj <- pracma::acot(g_adj)
  } else {
    alpha_adj <- atan(g_adj)
  }


  # Use previous identification before correction to obtain the right alpha
  if(names(alpha_hat) == "alpha_acot_ref" | names(alpha_hat) == "alpha_atan_ref") {
    alpha_adj <- pi - alpha_adj
  }

  list(
    "g_adj" = as.double(g_adj),
    "alpha_adj" = alpha_adj,
    "f_alpha" = f_alpha
  )

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

