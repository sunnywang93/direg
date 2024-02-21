
library(direg)
alpha_set <- c(0, pi/3, 5*pi/6)
N0 <- 150
M0 <- 101
sig <- 0.05
rout <- 320
H1 <- 0.8
H2 <- 0.5
delta_c <- 0.25
xout <- seq(0, 1, length.out = M0)
delta_grid <- seq(1/sqrt(M0), 0.4, length.out = 15)

Y_sum_learn <- purrr::map(seq_len(N0),
                         ~fbm_sheet(
                           t_n = M0,
                           e_n = M0,
                           alpha = alpha_set[2],
                           H1 = H1,
                           H2 = H2,
                           type = "sum",
                           sigma = sig
                           )
                         ) |>
  purrr::map(~list(t = xout,
                   X = .x))

image(
  Reduce('+', purrr::map(Y_sum_learn, ~.x$X^2)) / length(Y_sum_learn)
)


g_hat_sum <- estimate_angle(X_list = Y_sum_learn,
                            xout = xout,
                            delta = (1 / sqrt(M0)) * (1 + delta_c)
                            )

alpha_hat_sum <- identify_angle(angles = g_hat_sum,
                                X_list = Y_sum_learn,
                                dout = delta_grid,
                                xout = xout
                                )


u1 <- unname(c(cos(alpha_hat_sum), sin(alpha_hat_sum)))
u2 <- unname(c(-sin(alpha_hat_sum), cos(alpha_hat_sum)))

sigma_hat <- estimate_sigma(X_list = Y_sum_learn)

H_ani <- H_sheets_dir(X_list = Y_sum_learn,
                      tout = expand.grid(t1 = xout, t2 = xout),
                      delta = (1 / sqrt(M0)) * (1 + delta_c),
                      base_list = list(u1, u2),
                      sigma = sigma_hat)


delta_grid_L <- seq(0.01,
                    0.2,
                    length.out = 15)

test <- purrr::map_dbl(delta_grid_L,
                       ~L_sheets(X_list = Y_sum_learn,
                                 xout = xout,
                                 delta = .x,
                                 v = u1,
                                 Hv = H_ani[1]))

test2 <- purrr::map_dbl(delta_grid_L,
                        ~L_sheets(X_list = Y_sum_learn,
                                  xout = xout,
                                  delta = .x,
                                  v = u2,
                                  Hv = H_ani[2]))

L_sheets(X_list = Y_sum_learn,
         xout = xout,
         delta = 0.02,
         v = u2,
         Hv = H_ani[2])


L_sheets(X_list = Y_sum_learn,
         xout = xout,
         delta = (1 / sqrt(M0)) * (1 - 2*delta_c),
         v = u1,
         Hv = H_ani[1])






