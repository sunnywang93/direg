
#' Estimates the angle between two basis vectors
#'
#' @param X_list List, containing the following elements:
#' - **$t** Vector of sampling points,
#' - **$X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param xout Vector, containing the evaluation points at which
#' the angle should be computed.
#' @param delta Numeric, determining the spacings.
#' @export

estimate_angle <- function(X_list, xout, delta) {

  tout <- expand.grid(xout, xout)

  H_list <- H_sheets(X_list = X_list,
                     tout = tout,
                     delta = delta)


  tan_alpha <- mean(H_list$theta_e2 / H_list$theta_e1)**(1 / (2 * mean(H_list$H)))

  list(alpha_cot = pracma::acot(tan_alpha),
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



gaussian_blur <- function(matrix_data, sigma) {
  # Create a Gaussian kernel
  kernel_size <- ceiling(6 * sigma)
  if (kernel_size %% 2 == 0) kernel_size <- kernel_size + 1

  x <- seq(-3 * sigma, 3 * sigma, length.out = kernel_size)
  kernel <- dnorm(x, mean = 0, sd = sigma)
  kernel <- outer(kernel, kernel)
  kernel <- kernel / sum(kernel)

  # Apply the Gaussian kernel using convolution
  blurred_matrix <- matrix(0, nrow = nrow(matrix_data), ncol = ncol(matrix_data))

  for (i in 1:nrow(matrix_data)) {
    for (j in 1:ncol(matrix_data)) {
      min_row <- max(1, i - floor(kernel_size / 2))
      max_row <- min(nrow(matrix_data), i + floor(kernel_size / 2))
      min_col <- max(1, j - floor(kernel_size / 2))
      max_col <- min(ncol(matrix_data), j + floor(kernel_size / 2))

      neighborhood <- matrix_data[min_row:max_row, min_col:max_col]
      blurred_matrix[i, j] <- sum(neighborhood * kernel[1:(max_row - min_row + 1), 1:(max_col - min_col + 1)])
    }
  }

  return(blurred_matrix)
}

