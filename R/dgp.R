#' Simulate 1D fractional brownian motion
#'
#' Performs simulation of one-dimensional fractional
#' brownian motion on an equally spaced grid using fast-fourier transforms.
#'
#' @param H Numeric, the Hurst index of the fractional
#' brownian motion.
#' @param n Numeric, the number of evaluation points.
#' @param grid_max Numeric, the right-end point of the sampling interval.
#' The observed points are simulated on the interval `[0, grid_max]`, using
#' the self-similar property of fbms.
#' @return Vector, the observed points of the fbm.
#' @export


fbm_fft <- function(H, n, grid_max) {

  g <- ceiling(log(2*(n - 1)) / log(2))

  lambda <- -1

  while(sum(lambda < 0) >= 1) {
    m <- 2**g
    # Construct covariance vector - now based on n equally spaced points on
    # grid - how can we adapt the covariance to be evaluated at input grid points?
    cov_vec <- c()
    cov_vec[1] <- 1 # Variance of increments as the initial value
    for(k in 1:(m/2)) {
      cov_vec[k+1] <- 0.5 * ((k + 1)**(2*H) - 2*k**(2*H) + (k-1)**(2*H))
    }

    # Construct circulant vector
    cov_vec <- c(cov_vec, cov_vec[seq(length(cov_vec) - 1, 2)])

    # Calculate eigenvalues using fft
    lambda <- Re(fft(cov_vec))

    # Iterate g if the eigenvalues are not positive
    g <- g + 1

  }


  # Compute vector to perform fft on
  a <- c()
  a[1] <- sqrt(lambda[1] / m) * rnorm(1)
  a[2:(m/2 + 1)] <- sqrt(lambda[2:(m/2+1)] / (2 * n)) *
    (rnorm(m/2) + rnorm(m/2)*1i)
  a[m:(m/2 + 2)] <- sqrt(lambda[1:(m/2 - 1)] / (2 * n)) *
    (rnorm(m/2-1) - rnorm(m/2-1)*1i)

  # Obtain the differenced process
  W <- Re(fft(a))
  # Obtain the process at levels - self similarity
  # Scale by (T/n)**2H if we want to simulate on [0, T] instead of [0, 1]
  (grid_max / n)**(2 * H) * cumsum(W[1:n])

}


#' Simulates fractional brownian sheet
#'
#' Simulates a fractional brownian sheet on the two-dimensional canonical basis
#' based on the product of two fractional brownian motions (fbms) on another basis.
#' The individual fbms are simulated using fast fourier transforms for increased
#' speed; see `fbm_fft`.
#'

#' @param t1_n Numeric, containing the number of evaluation points for the
#' output along the direction of `u1`.
#' @param t2_n Numeric, containing the number of evaluation points for the
#' output along the direction of `u2`.
#' @param e1_n Numeric, containing the number of evaluation points along
#' the direction of the first canonical basis.
#' @param e2_n Numeric, containing the number of evaluation points along
#' the direction of the second canonical basis.
#' @param H1 Numeric, indicating the Hölder exponent of the first fractional
#' brownian motion.
#' @param H2 Numeric, indicating the Hölder exponent of the second fractional
#' brownian motion.
#' @returns Matrix, containing the fractional brownian sheets
#' observed on the canonical bases.
#' @export
fbm_sheet <- function(t_n, e_n, alpha, H1, H2) {

  # Specific coordinates of canonical bases
  e1 <- c(1, 0)
  e2 <- c(0, 1)

  # Specify coordinates of other bases
  u1 <- c(cos(alpha), sin(alpha))
  u2 <- c(-sin(alpha), cos(alpha))

  # Construct the evaluation points along the canonical bases
  t1_tilde <- seq(0, 1, length.out = e_n)
  t2_tilde <- seq(0, 1, length.out = e_n)

  # Construct the evaluation points along the other bases
  t1 <- seq(0 , cos(alpha) + sin(alpha), length.out = t_n)
  t2 <- seq(-sin(alpha), cos(alpha), length.out = t_n)

  # Simulate the 1D fractional brownian motions
  B1 <- fbm_fft(H = H1, n = length(t1), grid_max = cos(alpha) + sin(alpha))
  B2_tilde <- fbm_fft(H = H2, n = length(t1), grid_max = cos(alpha) + sin(alpha))

  # Obtain the indexes of the corresponding value of the fractional
  # brownian motion on the negative domain
  t2_idx_minus <- sapply(cos(alpha) - t2[t2 < 0],
                         function(x) which.min(abs(x - t1)))

  # Extract the values of the fbm on the negative part of the domain by
  # the self-similarity property
  B2_minus <- -B2_tilde[t2_idx_minus] +
    B2_tilde[which.min(abs( cos(alpha) - t1 ))]
  # Extract the values of the positive part
  B2_plus <- B2_tilde[seq_len(length(t1) - length(B2_minus))]
  # Combine them to get the full process
  B2 <- c(B2_minus, B2_plus)

  # Find the coordinates of the evaluation points defined on the
  # canonical basis with respect to the other basis
  t1_u1 <- apply(
    expand.grid(t2_tilde = t2_tilde, t1_tilde = t1_tilde),
    1,
    function(x) crossprod(x, c(cos(alpha), sin(alpha)))
  )

  t2_u2 <- apply(expand.grid(t2_tilde = t2_tilde, t1_tilde = t1_tilde),
                 1,
                 function(x) crossprod(x, c(-sin(alpha), cos(alpha)))
                 )

  # Obtain the indexes of the new coordinates
  t1_u1_idx <- sapply(t1_u1, function(x) which.min(abs(x - t1)))
  t2_u2_idx <- sapply(t2_u2, function(x) which.min(abs(x - t2)))

  # Extract the relevant values and build the 2D process
  X <- B1[t1_u1_idx] * B2[t2_u2_idx]

  # Convert the sequence into a matrix
  matrix(X,
         nrow = length(t1_tilde),
         ncol = length(t2_tilde))


}



fbm_prod <- function(H1, H2, n, endpoint) {

  fbm_1 <- fbm_fft(H1, n, endpoint)
  fbm_2 <- fbm_fft(H2, n, endpoint)

  apply(expand.grid(fbm_1, fbm_2), 1, function(x) x[1] * x[2]) |>
    matrix(nrow = n,
           ncol = n)

}




