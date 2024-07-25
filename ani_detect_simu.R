library(direg)
library(foreach)
library(doSNOW)
library(tictoc)
library(ggplot2)
library(here)


# Parameter settings
alpha_true <- pi/3 + pi/4
H_min <- 0.3
Hmax_set <- c(0.3, 0.4, 0.8)
N <- 100
M <- 101
sigma <- 0.05
delta_grid <- seq(1/sqrt(M), 0.4, length.out = 15)
rout <- 500



set.seed(123)
seeds <- sample.int(10000,
                    size = rout)

inter_dir <- here('intermediate')
result_folder <- here('result')
if(!dir.exists(result_folder)) {
  dir.create(result_folder)
}



if(!dir.exists(inter_dir)) {
  dir.create(inter_dir)
}

n_cores <- 50
# Create cluster nodes
cl <- makeCluster(spec = n_cores)
# Register cluster
registerDoSNOW(cl)

result <- foreach(H_max = Hmax_set, .combine = cbind) %:%
  foreach(i = 1:rout,
          .packages = c("direg"),
          .combine = 'c') %dopar%
    {

    each_filename <- paste0('result_',
                            as.character(i),
                            '.rda')

    each_filepath <- file.path(inter_dir,
                               each_filename)

    if (file.exists(each_filepath)) {
      next
    }

    set.seed(seeds[i])

    X_list_sum <- purrr::map(seq_len(N),
                             ~fbm_sheet(
                               t_n = M,
                               e_n = M,
                               alpha = alpha_true,
                               H1 = H_min,
                               H2 = H_max,
                               type = "sum",
                               sigma = sigma)
                             )
    # Check the variance of the simulated process
    #image(Reduce('+', purrr::map(X_list_sum, ~.x^2)) / length(X_list_sum))
    delta_grid <- seq(1/sqrt(M), 0.4, length.out = 10)

    xout <- seq(0, 1, length.out = M)
    sheets_list_sum <- purrr::map(X_list_sum,
                                  ~list(t = xout,
                                        X = .x))

    delta <- (1 / sqrt(M))

    # Estimation of angle
    alpha_sheet_sum <- estimate_angle(X_list = sheets_list_sum,
                                      xout = xout,
                                      delta = delta)

    #  Identification of angle
    alpha_unique_sum <- identify_angle(angles = alpha_sheet_sum,
                                       X_list = sheets_list_sum,
                                       dout = delta_grid,
                                       xout = xout)

    # Perform correction for remainder term of g
    alpha_hat_adj <- angle_correct(g_hat = alpha_sheet_sum$g_hat,
                                   alpha_hat = alpha_unique_sum$alpha,
                                   delta = delta,
                                   Hmax = alpha_unique_sum$H_max,
                                   Hmin = alpha_sheet_sum$H_min)

    # Perform anisotropic detection
    ani_indic <- ani_detect(X_list = sheets_list_sum,
                            beta_n = 10,
                            g_adj = alpha_hat_adj$g_adj,
                            alpha_adj = alpha_hat_adj$alpha_adj,
                            H_min = alpha_sheet_sum$H_min,
                            H_max = alpha_unique_sum$H_max,
                            delta = delta,
                            sigma = alpha_unique_sum$sigma)

    ani_indic

    }

colnames(result) <- paste0("H_", Hmax_set)

saveRDS(result,
        file = paste0(result_folder, "/H_max/result_all.rds")
        )

pct_0.3 <- sum(result_all[, "H_0.3"]) / rout * 100
pct_0.4 <- sum(result_all[, "H_0.4"]) / rout * 100
pct_0.8 <- sum(result_all[, "H_0.8"]) / rout * 100

