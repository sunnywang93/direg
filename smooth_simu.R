library(direg)
library(foreach)
library(parallel)
library(doSNOW)
library(ggplot2)
library(here)



alpha_set <- c(pi/30, pi/6, pi/4, pi/3, pi/2 - pi/30)
N0 <- 120
M0 <- 101
sigma_set <- c(0, 0.05, 0.1)
rout <- 200
H1 <- 0.8
H2 <- 0.5
xout <- seq(0, 1, length.out = M0)
delta_grid <- seq(1/sqrt(M0), 0.4, length.out = 15)
delta_learn <- (1 / sqrt(M0))


# Parameter settings for online set
M0_new_dense <- 151
M0_obs <- 101
xout_dense <- seq(0, 1, length.out = M0_new_dense)
xout_obs <- seq(0, 1, length.out = M0_obs)

# All configurations in loop
param_cart <- expand.grid(alpha = alpha_set,
                          sigma = sigma_set)


result_folder <- here('result_smooth')
if(!dir.exists(result_folder)) {
  dir.create(result_folder)
}

set.seed(1234)
seeds <- sample.int(10000, size = rout)

# Function for learning set generation and estimation
# For parallelization purposes

learn_simu <- function(k) {
  # Set varying parameters
  alpha <- param_cart[k, "alpha"]
  sigma <- param_cart[k, "sigma"]

  # Generate seeds for learning and online set
  set.seed(1234)

  # Generate learning set
  Y_sum_learn <- purrr::map(seq_len(N0),
                            ~fbm_sheet(
                              t_n = M0,
                              e_n = M0,
                              alpha = alpha,
                              H1 = H1,
                              H2 = H2,
                              type = "sum",
                              sigma = sigma
                            )
  ) |>
    purrr::map(~list(t = xout,
                     X = .x))

  # Estimate & identify angle from the learning set
  g_hat_sum <- estimate_angle(X_list = Y_sum_learn,
                              xout = xout,
                              delta = delta_learn)

  alpha_hat_sum <- identify_angle(angles = g_hat_sum,
                                  X_list = Y_sum_learn,
                                  dout = delta_grid,
                                  xout = xout)

  # Perform correction for estimated angle
  alpha_hat_sum_adj <- angle_correct(g_hat = g_hat_sum$g_hat,
                                     alpha_hat = alpha_hat_sum$alpha,
                                     delta = delta_learn,
                                     Hmax = alpha_hat_sum$H_max,
                                     Hmin = g_hat_sum$H_min)

  # Build anisotropic basis vectors
  u1 <- unname(c(cos(alpha_hat_sum_adj$alpha), sin(alpha_hat_sum_adj$alpha)))
  u2 <- unname(c(-sin(alpha_hat_sum_adj$alpha), cos(alpha_hat_sum_adj$alpha)))

  sigma_hat <- estimate_sigma(X_list = Y_sum_learn)

  H_ani <- H_sheets_dir(X_list = Y_sum_learn,
                        tout = expand.grid(t1 = xout, t2 = xout),
                        delta = delta_learn,
                        base_list = list(u1, u2),
                        sigma = sigma_hat)

  # Compute rotation matrix
  R_alpha <- matrix(c(cos(alpha_hat_sum_adj$alpha_adj),
                      -sin(alpha_hat_sum_adj$alpha_adj),
                      sin(alpha_hat_sum_adj$alpha_adj),
                      cos(alpha_hat_sum_adj$alpha_adj)
  ), nrow = 2, ncol = 2)


  # Comparison with KASSI 2023 - Estimate H and L along the canonical basis
  e1 <- c(1, 0)
  e2 <- c(0, 1)

  H_iso <- H_sheets_dir(X_list = Y_sum_learn,
                        tout = expand.grid(t1 = xout, t2 = xout),
                        delta = delta_learn,
                        base_list = list(e1, e2),
                        sigma = sigma_hat)

  list('R_alpha' = R_alpha,
       'H_ani' = H_ani,
       'H_iso' = H_iso,
       'alpha' = alpha,
       'sigma' = sigma)


}

# Set number of cores to use
n_cores <- 40
# Create cluster nodes
cl <- makeCluster(spec = n_cores)
# Register cluster
registerDoSNOW(cl)


result_list <- foreach(k = seq_len(nrow(param_cart))) %do% {
  # Run learning set operations and obtain results
  learn_list <- learn_simu(k = k)

  foreach(i = 1:rout,
          .packages = c("direg", "parallel", "snow"),
          .combine = 'c',
          .multicombine = TRUE)  %dopar% {

    set.seed(seeds[i])

    # Generate true curve on a dense grid, first without noise, then
    # manually add noise to the discretised ones

    # First without noise
    Y_new_true <- purrr::map(seq_len(1),
                             ~fbm_sheet(
                               t_n = M0_new_dense,
                               e_n = M0_new_dense,
                               alpha = learn_list$alpha,
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
                                                   sd = learn_list$sigma)
                              )
    )

    # Now construct the smoothing grid - we need to transform it onto the
    # anisotropic basis!
    tout <- expand.grid(t1 = xout_dense,
                        t2 = xout_dense)

    # Rotation on the smoothing grid
    tout_rot <- t(apply(tout, 1,
                        function(x) crossprod(t(learn_list$R_alpha), x)))

    # Perform rotation on the observed grid
    tobs_rot <- t(apply(tobs, 1,
                        function(x) crossprod(t(learn_list$R_alpha), x)))

    # Compute optimal smoothing bandwidth
    h_star <- bw_smooth(H1 = learn_list$H_ani[1],
                        H2 = learn_list$H_ani[2],
                        M0 = M0_obs^2)

    # Compute smoothing bandwidth along canonical basis (w/o directional reg)
    h_star_iso <- bw_smooth(H1 = learn_list$H_iso[1],
                            H2 = learn_list$H_iso[2],
                            M0 = M0_obs^2)

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

    risk_rel <- risk_ani / risk_iso

    list(
      list('risk_ani' = risk_ani,
           'risk_iso' = risk_iso,
           'risk_rel' = risk_rel,
           'h1_ani' = h_star[1],
           'h2_ani' = h_star[2],
           'h1_iso' = h_star_iso[1],
           'h2_iso' = h_star_iso[2])
    )

  }

}

saveRDS(result_list,
        file = paste0(here(), "/result_smooth/result_list.rds"))


# analysis of results =========================================================
library(dplyr)
library(tidyr)
library(stringr)

list_names <- purrr::map_chr(seq_len(nrow(param_cart)),
                             ~paste0("alpha",
                                     round(param_cart[.x, "alpha"], 2),
                                     "_sigma",
                                     param_cart[.x, "sigma"])
                             )

names(result_list) <- list_names

result_df <- bind_rows(
  lapply(seq_along(result_list), function(result_id) {
    df_result <- as.data.frame(do.call(rbind, result_list[[result_id]]))
    as.data.frame(
      cbind(alpha = as.numeric(str_extract(list_names[result_id],
                                           "(?<=alpha)[0-9.]+")),
            sigma = as.numeric(str_extract(list_names[result_id],
                                           "(?<=sigma)[0-9.]+")),
            risk_ani = unlist(df_result$risk_ani),
            risk_iso = unlist(df_result$risk_iso),
            risk_rel = unlist(df_result$risk_rel),
            h1_ani = purrr::map_dbl(df_result$h_ani, ~.x[1]),
            h2_ani = purrr::map_dbl(df_result$h_ani, ~.x[2]),
            h1_iso = purrr::map_dbl(df_result$h_iso, ~.x[1]),
            h2_iso = purrr::map_dbl(df_result$h_iso, ~.x[2]))
    )
  })
)

result_df |>
  select(alpha, sigma, risk_rel) |>
  ggplot(
    aes(x = as.factor(sigma), y = risk_rel, fill = as.factor(alpha))
  ) +
  geom_boxplot() +
  geom_hline(yintercept = 1, col = "red")


result_df |>
  mutate(log_risk = log(risk_ani / risk_iso)) |>
  ggplot(aes(x = as.factor(sigma), y = log_risk, fill = as.factor(alpha))) +
  geom_boxplot()


