
library(direg)

# Parameter settings
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


# Generate learning set
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

# Check variance
image(
  Reduce('+', purrr::map(Y_sum_learn, ~.x$X^2)) / length(Y_sum_learn)
)

# Estimate & identify angle from the learning set
g_hat_sum <- estimate_angle(X_list = Y_sum_learn,
                            xout = xout,
                            delta = (1 / sqrt(M0)) * (1 + delta_c)
                            )

alpha_hat_sum <- identify_angle(angles = g_hat_sum,
                                X_list = Y_sum_learn,
                                dout = delta_grid,
                                xout = xout
                                )

# Estimate other parameters required for smoothing

# Build anisotropic basis vectors
u1 <- unname(c(cos(alpha_hat_sum), sin(alpha_hat_sum)))
u2 <- unname(c(-sin(alpha_hat_sum), cos(alpha_hat_sum)))

sigma_hat <- estimate_sigma(X_list = Y_sum_learn)

H_ani <- H_sheets_dir(X_list = Y_sum_learn,
                      tout = expand.grid(t1 = xout, t2 = xout),
                      delta = (1 / sqrt(M0)) * (1 + delta_c),
                      base_list = list(u1, u2),
                      sigma = sigma_hat)

L_ani <- c(L_sheets(X_list = Y_sum_learn,
                    xout = xout,
                    delta = (1 / sqrt(M0)) * (1 - 2*delta_c),
                    v = u1,
                    Hv = H_ani[1]),
           L_sheets(X_list = Y_sum_learn,
                    xout = xout,
                    delta = (1 / sqrt(M0)) * (1 - 2*delta_c),
                    v = u2,
                    Hv = H_ani[2])
           )


# Compute rotation matrix
R_alpha <- matrix(c(cos(alpha_hat_sum), -sin(alpha_hat_sum),
                    sin(alpha_hat_sum), cos(alpha_hat_sum)
                    ), nrow = 2, ncol = 2)


# Comparison with KASSI 2023 - Estimate H and L along the canonical basis
e1 <- c(1, 0)
e2 <- c(0, 1)

H_iso <- H_sheets_dir(X_list = Y_sum_learn,
                      tout = expand.grid(t1 = xout, t2 = xout),
                      delta = (1 / sqrt(M0)) * (1 + delta_c),
                      base_list = list(e1, e2),
                      sigma = sigma_hat)

L_iso <- c(L_sheets(X_list = Y_sum_learn,
                    xout = xout,
                    delta = (1 / sqrt(M0)) * (1 - 2*delta_c),
                    v = e1,
                    Hv = H_iso[1]),
           L_sheets(X_list = Y_sum_learn,
                    xout = xout,
                    delta = (1 / sqrt(M0)) * (1 - 2*delta_c),
                    v = e2,
                    Hv = H_iso[2])
)


# ============================================== Finish learning set


# =========================================== Testing set
M0_new_dense <- 201
M0_obs <- 101
xout_dense <- seq(0, 1, length.out = M0_new_dense)
xout_obs <- seq(0, 1, length.out = M0_obs)
k <- (3/ (4 * 0.99))^2

# Generate true curve on a dense grid, first without noise, then
# manually add noise to the discretised ones
Y_new_true <- purrr::map(seq_len(20),
                         ~fbm_sheet(
                           t_n = M0_new_dense,
                           e_n = M0_new_dense,
                           alpha = alpha_set[2],
                           H1 = H1,
                           H2 = H2,
                           type = "sum",
                           sigma = 0)
                         ) |>
  purrr::map(~list(t = xout_dense,
                   X = .x))

# Discretise the true surface with interpolation
tobs <- expand.grid(t1 = xout_obs,
                    t2 = xout_obs)

Y_new_obs <- purrr::map(Y_new_true,
                        ~pracma::interp2(x = .x$t,
                                         y = .x$t,
                                         Z = .x$X,
                                         xp = tobs[, 2],
                                         yp = tobs[, 1],
                                         method = "nearest") |>
                          matrix(nrow = M0_obs,
                                 ncol = M0_obs)
                        )

# Now add some noise
Y_new_noisy <- purrr::map(Y_new_obs,
                          ~list(t = xout_obs,
                                X = .x + rnorm(n = M0_obs^2,
                                               sd = sig)
                                )
                          )

# Now construct the smoothing grid - we need to transform it onto the
# anisotropic basis!
tout <- expand.grid(t1 = xout_dense,
                    t2 = xout_dense)

tout_rot <- t(apply(tout, 1, function(x) crossprod(t(R_alpha), x)))

# Perform rotation on the observed grid
tobs_rot <- t(apply(tobs, 1, function(x) crossprod(t(R_alpha), x)))

# Compute optimal smoothing bandwidth
h_star <- bw_smooth(H1 = H_ani[1],
                    H2 = H_ani[2],
                    L1 = L_ani[1],
                    L2 = L_ani[2],
                    sigma = sigma_hat,
                    k = k,
                    M0 = M0_obs^2,
                    rate = TRUE)

# Compute smoothing bandwidth along canonical basis (w/o directional reg)
h_star_iso <- bw_smooth(H1 = H_iso[1],
                        H2 = H_iso[2],
                        L1 = L_iso[1],
                        L2 = L_iso[2],
                        sigma = sigma_hat,
                        k = k,
                        M0 = M0_obs^2,
                        rate = TRUE)

# Smooth surfaces
Y_smoothed <- ksmooth_bi(Y_list = Y_new_noisy,
                         bw_vec = h_star,
                         tobs = tobs_rot,
                         tout = tout_rot)

Y_smoothed_iso <- ksmooth_bi(Y_list = Y_new_noisy,
                             bw_vec = h_star_iso,
                             tobs = tobs,
                             tout = tout)

# Compute the risk
risk_ani <- purrr::map2_dbl(Y_new_true, Y_smoothed,
                            ~mean(abs(.x$X - .y)))


risk_iso <- purrr::map2_dbl(Y_new_true, Y_smoothed_iso,
                            ~mean(abs(.x$X - .y)))


par(mfrow = c(1, 4))
image(Y_new_true[[11]]$X)
image(Y_new_noisy[[11]]$X)
image(Y_smoothed[[11]])
image(Y_smoothed_iso[[11]])



