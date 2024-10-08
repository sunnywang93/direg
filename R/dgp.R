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
    m <- 2**g - 1
    # Construct covariance vector - now based on n equally spaced points on
    # grid - how can we adapt the covariance to be evaluated at input grid points?
    cov_vec <- c()
    cov_vec[1] <- 1 # Variance of increments as the initial value
    for(k in 1:((m+1) / 2)) {
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
  a[(m+1)/2 + 1] <- sqrt(lambda[(m+1)/2 + 1] / m) * rnorm(1)

  U <- rnorm(n = (m+1)/2 - 1)
  V <- rnorm(n = (m+1)/2 - 1)

  a[2:((m+1) / 2)] <- sqrt(lambda[2:((m+1) / 2)] / (2 * m)) *
    (U + V*1i)

  a[(m+1):((m+1)/2 + 2)] <- sqrt(lambda[2:((m+1) / 2)] / (2 * m)) *
    (U - V*1i)

  # Obtain the differenced process
  W <- Re(fft(a))
  # Obtain the process at levels - self similarity
  # Scale by (T/n)**2H if we want to simulate on [0, T] instead of [0, 1]
  (grid_max / n)**(1 * H) * cumsum(W[1:n])

}

#' Simulates fractional brownian sheet with constant anisotropy
#'
#' Simulates a fractional brownian sheet on the two-dimensional canonical basis
#' based on the product of two fractional brownian motions (fbms) on another basis.
#' In order to avoid confusions with interpretation, `H1` should always be greater
#' than `H2`. The individual fbms are simulated using fast fourier transforms for increased
#' speed; see `fbm_fft`.
#'

#' @param t_n Numeric, the number of evaluation points for the
#' output along the direction of the non-canonical bases.
#' @param e_n Numeric, the number of evaluation points along
#' the direction of the canonical bases
#' @param alpha Numeric, the angle between the first non-canonical and canonical
#' basis. Must be a number between `[0, pi/2]`.
#' @param H1 Numeric, the Hölder exponent of the first fractional
#' brownian motion.
#' @param H2 Numeric, the Hölder exponent of the second fractional
#' brownian motion.
#' @param type String, either "sum" or "prod", indicating whether to simulate
#' from the sum or products of two fBms.
#' @returns Matrix, containing the fractional brownian sheets
#' observed on the canonical bases.
#' @export
fbm_sheet <- function(t_n, e_n, alpha, H1, H2, type = "sum", sigma = NULL) {

  # Specific coordinates of canonical bases
  e1 <- c(1, 0)
  e2 <- c(0, 1)

  if(H1 <= H2) {
    alpha <- alpha + pi/2
  }

  while(alpha >= pi) {
    alpha <- alpha - pi
  }

  # if(alpha >= pi) {
  #   alpha <- alpha - pi
  # }

  # Specify coordinates of other bases
  u1 <- c(cos(alpha), sin(alpha))
  u2 <- c(-sin(alpha), cos(alpha))

  # Construct the evaluation points along the canonical bases
  t1_tilde <- seq(0, 1, length.out = e_n)
  t2_tilde <- seq(0, 1, length.out = e_n)

  # Construct the evaluation points along the other bases
  t1 <- seq(0, abs(cos(alpha)) + sin(alpha), length.out = t_n)
  # Simulate the 1D fractional brownian motions
  if(alpha <= pi/2) {
    t2 <- seq(-sin(alpha), cos(alpha), length.out = t_n)
    B1 <- fbm_fft(H = H1, n = length(t1), grid_max = max(t1))

    B2_tilde <- fbm_fft(H = H2, n = length(t1),
                        grid_max = max(t1))
    # Obtain the indexes of the corresponding value of the fractional
    # brownian motion on the negative domain
    t2_idx_minus <- sapply(-t2[t2 < 0], function(x) which.min(abs(x - t1)))

    # Extract the values of the fbm on the negative part of the domain by
    # the self-similarity property
    if(length(t2_idx_minus) == 0) {
      B2_minus <- c()
    } else {
      B2_minus <- -B2_tilde[t2_idx_minus]
    }
    # Extract the values of the positive part
    B2_plus <- B2_tilde[seq_len(length(t1) - length(B2_minus))]
    # Combine them to get the full process
    B2 <- c(B2_minus, B2_plus)

  } else {
    t1_proj <- seq(cos(alpha), sin(alpha), length.out = t_n)
    t2 <- seq(cos(alpha) - sin(alpha), 0, length.out = t_n)

    B1_tilde <- fbm_fft(H = H1, n = length(t1), grid_max = max(t1))

    t1_idx_minus <- sapply(-t1_proj[t1_proj < 0],
                           function(x) which.min(abs(x - t1)))

    B1_minus <- -B1_tilde[t1_idx_minus]

    B1_plus <- B1_tilde[seq_len(length(t1) - length(B1_minus))]
    B1 <- c(B1_minus, B1_plus)

    B2_tilde <- fbm_fft(H = H2, n = length(t1), grid_max = max(t1))

    t2_idx_minus <- sapply(-t2[t2 < 0],
                           function(x) which.min(abs(x - t1)))

    B2_minus <- -B2_tilde[t2_idx_minus]
    # Extract the values of the positive part
    B2_plus <- B2_tilde[seq_len(length(t2) - length(B2_minus))]
    # Combine them to get the full process
    B2 <- c(B2_minus, B2_plus)


  }

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
  if(type == "sum") {
    X <- B1[t1_u1_idx] + B2[t2_u2_idx]
  } else {
    X <- B1[t1_u1_idx] * B2[t2_u2_idx]
  }


  # Convert the sequence into a matrix
  if(missing(sigma)) {
    matrix(X,
           nrow = length(t1_tilde),
           ncol = length(t2_tilde))
  } else {
    noise <- rnorm(n = length(X), mean = 0, sd = sigma)
    matrix(X + noise,
           nrow = length(t1_tilde),
           ncol = length(t2_tilde))
  }

}


#' Simulates fractional brownian sheet with time-varying anisotropy
#'
#' Simulates a fractional brownian sheet on the two-dimensional canonical basis
#' based on the sum of two fractional brownian motions (fbms) on another basis.
#' The individual fbms are simulated using fast fourier transforms for increased
#' speed; see `fbm_fft`. This function simulates brownian sheets with time-varying
#' angles; for a constant angle, please use `fbm_sheet` instead.
#'

#' @param t_n Numeric, containing the number of evaluation points for the
#' output along the direction of anisotropic bases.
#' @param e_n Numeric, containing the number of evaluation points along
#' the direction of the canonical bases.
#' @param alpha_fun Function, indicating how to construct the time-varying angles.
#' @param H1 Numeric, indicating the Hölder exponent of the first fractional
#' brownian motion.
#' @param H2 Numeric, indicating the Hölder exponent of the second fractional
#' brownian motion.
#' @param b_n Numeric, indicating the number of evaluation points to simulate
#' the fBms on.
#' @param type String, either "sum" or "prod", indicating whether to simulate
#' from the sum or products of two fBms.
#' @returns Matrix, containing the fractional brownian sheets
#' observed on the canonical bases.
#' @export
fbm_sheet_var <- function(t_n, e_n, alpha_fun, H1, H2, b_n = 10000,
                          type = "sum") {

  # Specific coordinates of canonical bases
  e1 <- c(1, 0)
  e2 <- c(0, 1)

  # Construct the evaluation points along the canonical bases
  t1_tilde <- seq(0, 1, length.out = e_n)
  t2_tilde <- seq(0, 1, length.out = e_n)

  # Construct the 2D grid
  tt <- expand.grid(t1 = t1_tilde, t2 = t2_tilde)

  # Construct alpha(t)
  #alpha_t <- alpha_fun(tt$t1, tt$t2)
  if(length(formals(alpha_fun)) == 2) {
    alpha_t <- mapply(alpha_fun, tt$t1, tt$t2)
  } else {
    alpha_t <- apply(tt, 1, alpha_fun)
  }

  t1 <- sapply(seq_along(alpha_t),
               function(idx) (tt$t1[idx] * cos(alpha_t[idx])) +
                 tt$t2[idx] * sin(alpha_t[idx]))


  t2 <- sapply(seq_along(alpha_t),
               function(idx) (-tt$t1[idx] * sin(alpha_t[idx])) +
                 tt$t2[idx] * cos(alpha_t[idx]))

  # Simulate the 1D fractional brownian motions
  B1_tilde <- fbm_fft(H = H1, n = b_n, grid_max = max(t1))

  t1_grid_max <- seq(0, max(t1), length.out = b_n)

  t2_grid_max <- seq(0, max(t2) - min(t2), length.out = b_n)

  B2_tilde <- fbm_fft(H = H2, n = b_n, grid_max = max(t2) - min(t2))

  # Obtain the indexes of the corresponding values of the fractional
  # brownian motion on the negative domain
  t2_idx_minus <- sapply(max(t2) - t2[t2 < 0],
                         function(x) which.min(abs(x - t2_grid_max)))

  # Extract the values of the fbm on the negative part of the domain by
  # the self-similarity property
  B2_minus <- -B2_tilde[t2_idx_minus] +
    B2_tilde[which.min(abs(max(t2) - t2_grid_max))]

  # Extract the values of the positive part
  t2_idx_plus <- sapply(t2[t2 >= 0],
                        function(x) which.min(abs(x - t2_grid_max)))

  B2_plus <- B2_tilde[t2_idx_plus]
  # Combine them to get the full process

  # Find indexes of negative t2s in the original t2 vectors
  t2_minus_orig_idx <- which(t2 < 0)

  # Find indexes of positive t2s in the original t2 vectors
  t2_plus_orig_idx <- which(t2 >= 0)

  B2_named <- c(setNames(B2_minus, t2_minus_orig_idx),
                setNames(B2_plus, t2_plus_orig_idx))

  B2 <- B2_named[order(as.double(names(B2_named)))]

  t1_idx <- sapply(t1,
                   function(x) which.min(abs(x - t1_grid_max)))

  B1 <- B1_tilde[t1_idx]

  # Extract the relevant values and build the 2D process
  if(type == "sum") {
    X <- B1 + B2
  } else {
    X <- B1 * B2
  }

  # Convert the sequence into a matrix
  matrix(unname(X),
         nrow = length(t1_tilde),
         ncol = length(t2_tilde))

}


#' Simulate isotropic process based on fractional brownian motions

fbm_sum <- function(H1, H2, n, endpoint) {

  fbm_1 <- fbm_fft(H1, n, endpoint)
  fbm_2 <- fbm_fft(H2, n, endpoint)

  apply(expand.grid(fbm_1, fbm_2), 1, function(x) x[1] + x[2]) |>
    matrix(nrow = n,
           ncol = n)

}



