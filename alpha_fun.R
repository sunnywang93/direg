library(direg)
library(foreach)
library(parallel)
library(tictoc)
library(snow)
library(ggplot2)
library(here)

# Parameter settings
N <- 100
M <- 101
H1 <- 0.8
H2 <- 0.5
delta_grid <- seq(0.05, 0.4, length.out = 21)
rout <- 20
xtrue <- seq(0, 1, length.out = M)
xparam <- seq(0, 1, length.out = 21)


# Different functions for alpha

alpha_piece <- function(t1, t2) {
  if (t2 < 0.5) {
    if (t1 < 0.5) {
      return(pi / 6)
    } else {
      return(pi / 5)
    }
  } else {
    if (t1 < 0.5) {
      return(pi / 4)
    } else {
      return(pi / 3)
    }
  }
}


alpha_slow <- function(t) {
  log(sqrt(crossprod(t, t)) + 1) + (pi / 9)
}

alpha_piece_smooth <- function(t1, t2, npart = 5) {

  pi / (8 * npart) * (2*ceiling(npart * t1) + 2*ceiling(npart * t2) - 2)
}



X_list <- purrr::map(seq_len(N),
                     ~fbm_sheet_var(
                       t_n = M,
                       e_n = M,
                       alpha_fun = alpha_slow,
                       H1 = .8,
                       H2 = .5)
)

image(Reduce('+', purrr::map(X_list, ~.x^2)) / length(X_list))

image(matrix(apply(tt, 1, alpha_slow), nrow = M, ncol = M))


sheets_list <- purrr::map(X_list,
                          ~list(t = seq(0, 1, length.out = M),
                                X = .x))


alpha_sheet <- lapply(delta_grid, function(delta) {
  estimate_angle(X_list = sheets_list,
                 xout = xparam,
                 delta = delta)
})


alpha_avg_cot <- purrr::map(alpha_sheet,
                            ~alpha_avg(X = .x$alpha_cot,
                                       npart = 20))



# Function to plot a matrix as an image
plot_matrix_as_image <- function(matrix_data, title) {
  image(matrix_data, main = title)
}


# Iterate through the list of matrices and plot each one
for (i in 1:length(delta_grid)) {
  plot_matrix_as_image(alpha_sheet[[i]]$alpha_cot, paste("Delta", i))
}

for (i in 1:length(delta_grid)) {
  plot_matrix_as_image(alpha_avg_cot[[i]], paste("Delta", i))
}

sigma = 0.5
blurred_matrix <- gaussian_blur(alpha_sheet[[2]]$alpha_cot, sigma)


