library(direg)
library(foreach)
library(parallel)
library(tictoc)
library(snow)
library(ggplot2)
library(here)

# Parameter settings
N <- c(100, 150)
Mset <- c(51, 101)
H1 <- 0.8
H2 <- 0.5
rout <- 350
alpha_set <- c(pi/30, pi/6, pi/5, pi/4, pi/3, pi/2 - pi/30, 2*pi/3, 3*pi/4, 5*pi/6)
delta_c <- 0.25
sigma_set <- c(0, 0.025, 0.05, 0.1)


param_cart <- expand.grid(M = Mset,
                          alpha = alpha_set,
                          sigma = sigma_set)
# Set seeds to ensure reproducibility
set.seed(123)

seeds <- sample.int(10000,
                    size = rout)

inter_dir <- here('intermediate')
result_folder <- here('result')
if(!dir.exists(result_folder)) {
  dir.create(result_folder)
}


for(k in 1:nrow(param_cart)) {


  if(!dir.exists(inter_dir)) {
    dir.create(inter_dir)
  }
  # Set number of cores to use
  n_cores <- 70
  # Create cluster nodes
  cl <- makeCluster(spec = n_cores)
  # Register cluster
  doSNOW::registerDoSNOW(cl)
  # Create progress bar
  pb <- txtProgressBar(min = 1,
                       max = rout,
                       style = 3)

  opts <- list(progress = function(n) setTxtProgressBar(pb, n))

  alpha_true <- param_cart[k, "alpha"]
  M <- param_cart[k, "M"]
  delta_grid <- seq(1/sqrt(M), 0.4, length.out = 15)
  sigma <- param_cart[k, "sigma"]


  tic()
  foreach(i = 1:rout,
          .packages = c("direg"),
          .options.snow = opts) %dopar%
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
      # Generate sum of fBms
      X_list_sum <- purrr::map(seq_len(N),
                               ~fbm_sheet(
                                 t_n = M,
                                 e_n = M,
                                 alpha = alpha_true,
                                 H1 = H1,
                                 H2 = H2,
                                 type = "sum",
                                 sigma = sigma)
      )
      # Check the variance of the simulated process
      #image(Reduce('+', purrr::map(X_list_sum, ~.x^2)) / length(X_list_sum))
      delta_grid <- seq(1/sqrt(M), 0.4, length.out = 15)

      xout <- seq(0, 1, length.out = M)
      sheets_list_sum <- purrr::map(X_list_sum,
                                    ~list(t = xout,
                                          X = .x))

      delta <- (1 / sqrt(M)) * (1 + delta_c)

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

      # Compute true g for comparison purposes
      g_true <- g_compute(alpha = alpha_true,
                          H1 = H1,
                          H2 = H2,
                          delta = delta)


      # Perform iteration to ensure we compare to the right alpha (recall that
      # true alpha is defined up to k*pi/2 additions)
      while(alpha_true > pi) {
        alpha_true <- alpha_true - pi
      }


      result <- list(
        'alpha_hat' = alpha_unique_sum$alpha,
        'alpha_hat_adj' = alpha_hat_adj$alpha_adj,
        'ghat' = alpha_sheet_sum$g_hat,
        'ghat_adj' = alpha_hat_adj$g_adj,
        'g_true' = g_true,
        'f_alpha' = alpha_hat_adj$f_alpha
      )

      save(result,
           file = each_filepath)

    }

  stopCluster(cl)
  toc()

  fls <- list.files(inter_dir,
                    pattern = ".rda")


  result_list <- lapply(fls,
                        function(x) get(eval(load(paste0(inter_dir, '/', x
                        )))))


  saveRDS(result_list,
          file = paste0(result_folder, "/N", N,
                        "_M", M,
                        "_alpha", round(alpha_true, 2),
                        "_sigma", sigma,
                        ".rds")

  )

  fs::file_delete(inter_dir)

}

library(here)
folder_path <- here("result")
file_names <- list.files(path = folder_path, pattern = "\\.rds$",
                         full.names = TRUE)

dat_list <- lapply(file_names, function(x) readRDS(x)) |>
  lapply(function(dat) t(sapply(dat, function(x) x)))

names(dat_list) <- basename(file_names)

library(dplyr)
library(tidyr)
result_df <- bind_rows(
  lapply(names(dat_list), function(name) {
    mat <- dat_list[[name]]
    name_parts <- unlist(strsplit(name, "_"))
    data.frame(
      N = as.integer(sub("N(\\d+).*", "\\1", name_parts[1])),
      M = as.integer(sub("M(\\d+).*", "\\1", name_parts[2])),
      alpha = as.numeric(sub("alpha([0-9.]+).*", "\\1", name_parts[3])),
      sigma = as.numeric(sub("sigma([0-9.]+)\\.rds", "\\1", name_parts[4])),
      mat
    )
  })
)

for(i in 5:length(result_df)) {

  result_df[[i]] <- unlist(result_df[[i]])

}

result_long <- result_df |>
  pivot_longer(cols = c(alpha_adj_sum, risk_sum),
               names_to = c(".value", "dgp"),
               names_pattern = "(.+)_(.+)",
               values_to = "value")

# result_long <- result_long |>
#   mutate(risk = risk - pi/2)
#
#
# result_long <- result_long |>
#   mutate(risk = abs(risk))


boxplot_N100_M51_sum <- result_df |>
  filter(N == 100, M == 51) |>
  ggplot(aes(x = as.factor(sigma), y = risk_sum,
             fill = as.factor(alpha))) +
  geom_boxplot() +
  ggtitle("Boxplots $f_1 = B_1 + B2 (N = 100, M = 51)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\alpha}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1) +
  theme(plot.title = element_text(hjust = 0.5))


boxplot_N100_M101_sum <- result_df |>
  filter(N == 100, M == 101) |>
  ggplot(aes(x = as.factor(sigma), y = risk_sum,
             fill = as.factor(alpha))) +
  geom_boxplot() +
  ggtitle("Boxplots $f_1 = B_1 + B2 (N = 100, M = 101)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\alpha}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1) +
  theme(plot.title = element_text(hjust = 0.5))




test = result_df |>
  filter(N == 100, M == 101) |>
  mutate(risk_rel = risk_sum / abs(alpha_prelim_sum - alpha))

lim_test <- filter(test,
                   N == 100,
                   M == 101) |>
  pull(risk_rel) |>
  quantile(c(.05, .95))

test |>
  filter(N == 100, M == 101) |>
  mutate(risk_rel = risk_sum / abs(alpha_prelim_sum - alpha)) |>
  ggplot(aes(x = as.factor(sigma), y = risk_rel,
             fill = as.factor(alpha))) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = lim_test) +
  ggtitle("Boxplots $f_1 = B_1 + B2 (N = 100, M = 101)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\alpha}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1) +
  theme(plot.title = element_text(hjust = 0.5))


test |>
  filter(N == 100, M == 101) |>
  ggplot(aes(x = as.factor(sigma), y = abs(alpha_prelim_sum - alpha),
             fill = as.factor(alpha))) +
  geom_boxplot() +
  ggtitle("Boxplots $f_1 = B_1 + B2 (N = 100, M = 101)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\alpha}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1) +
  theme(plot.title = element_text(hjust = 0.5))


boxplot_N100_M101_sum_ghat <- result_df |>
  filter(N == 100, M == 101) |>
  ggplot(aes(x = as.factor(sigma), y = ghat_risk,
             fill = as.factor(alpha))) +
  geom_boxplot() +
  ggtitle("Boxplots $f_1 = B_1 + B2 (N = 100, M = 101)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\widehat g_{\\alpha}}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1) +
  theme(plot.title = element_text(hjust = 0.5))


boxplot_N100_M51_sum_ghat <- result_df |>
  filter(N == 100, M == 51) |>
  ggplot(aes(x = as.factor(sigma), y = ghat_risk,
             fill = as.factor(alpha))) +
  geom_boxplot() +
  ggtitle("Boxplots $f_1 = B_1 + B2 (N = 100, M = 51)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\widehat g_{\\alpha}}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1) +
  theme(plot.title = element_text(hjust = 0.5))

library(tikzDevice)

tikz(file = "/home/swang/directional_regularity/plots/alpha_explained2.tex",
     width = 5,
     height = 5,
     standAlone = TRUE)
plot(boxplot_N100_M51_sum)
dev.off()


lim_N100_M51_prod <- filter(result_long,
                            dgp == "prod",
                            N == 100,
                            M == 51) |>
  pull(risk) |>
  quantile(c(.1, .9))

boxplot_N100_M51_prod <- result_long |>
  ggplot(aes(x = as.factor(sigma), y = risk,
             fill = as.factor(alpha))) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = quantile(lim_N100_M51_prod, c(0.1, 0.9))) +
  ggtitle("Boxplots $f_1 = B_1 * B2 (N = 100, M = 51)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\alpha}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1)


boxplot_N100_M101_sum <- result_long |>
  filter(dgp == "sum", N == 100, M == 101) |>
  ggplot(aes(x = as.factor(sigma), y = risk,
             fill = as.factor(alpha))) +
  geom_boxplot() +
  ggtitle("Boxplots $f_1 = B_1 + B2 (N = 100, M = 101)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\alpha}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1)


lim_N100_N51_prod <- filter(result_long,
                            dgp == "prod",
                            N == 100,
                            M == 51) |>
  pull(risk) |>
  quantile(c(.15, .85))


boxplot_N100_M51_prod <- result_long |>
  filter(dgp == "prod", N == 100, M == 51) |>
  ggplot(aes(x = as.factor(sigma), y = risk,
             fill = as.factor(alpha))) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = lim_N100_N51_prod) +
  ggtitle("Boxplots $f_1 = B_1 * B2 (N = 100, M = 51)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\alpha}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1)


boxplot_N200_M26_sum <- result_long |>
  filter(dgp == "sum", N == 200, M == 26) |>
  ggplot(aes(x = as.factor(sigma), y = risk,
             fill = as.factor(alpha))) +
  geom_boxplot() +
  ggtitle("Boxplots $f_1 = B_1 + B2 (N = 200, M = 26)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\alpha}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1)


boxplot_N200_M51_sum <- result_long |>
  filter(dgp == "sum", N == 200, M == 51) |>
  ggplot(aes(x = as.factor(sigma), y = risk,
             fill = as.factor(alpha))) +
  geom_boxplot() +
  ggtitle("Boxplots $f_1 = B_1 + B2 (N = 200, M = 51)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\alpha}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1)


lim_N200_N51_prod <- filter(result_long,
                            dgp == "prod",
                            N == 200,
                            M == 51) |>
  pull(risk) |>
  quantile(c(.15, .85))

boxplot_N200_M51_prod <- result_long |>
  filter(dgp == "prod", N == 200, M == 51) |>
  ggplot(aes(x = as.factor(sigma), y = risk,
             fill = as.factor(alpha))) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = lim_N200_N51_prod) +
  ggtitle("Boxplots $f_1 = B_1 * B2 (N = 200, M = 51)$") +
  xlab("$\\sigma$") +
  ylab("$\\mathcal{R}_{\\alpha}$") +
  labs(fill = "$\\alpha$") +
  theme_minimal() +
  scale_fill_grey(start = 0.3, end = 1)


# library(tikzDevice)
# tikz(file = "N100_M26_sum.tex", width = 5, height = 5,
#      standAlone = TRUE)
# plot(boxplot_N100_M26_sum)
# dev.off()
#
#
# tikz(file = "N100_M51_sum.tex", width = 5, height = 5,
#      standAlone = TRUE)
# plot(boxplot_N100_M51_sum)
# dev.off()
#
# tikz(file = "N200_M26_sum.tex", width = 5, height = 5,
#      standAlone = TRUE)
# plot(boxplot_N200_M26_sum)
# dev.off()
#
# tikz(file = "N200_M51_sum.tex", width = 5, height = 5,
#      standAlone = TRUE)
# plot(boxplot_N200_M51_sum)
# dev.off()
#
#
# tikz(file = "N100_M51_prod.tex", width = 5, height = 5,
#      standAlone = TRUE)
# plot(boxplot_N100_M51_prod)
# dev.off()
#
# tikz(file = "N200_M51_prod.tex", width = 5, height = 5,
#      standAlone = TRUE)
# plot(boxplot_N200_M51_prod)
# dev.off()
#
#
# tikz(file = "var_pi4.tex", width = 5, height = 5,
#      standAlone = TRUE)
#
# image(Reduce('+', purrr::map(X_list_prod, ~.x^2)) /
#         length(X_list_prod),
#       main = "Variance of Anisotropic process",
#       xlab = "t1",
#       ylab = "t2")
# dev.off()
#
#
# tikz(file = "iso.tex", width = 5, height = 5,
#      standAlone = TRUE)
# image(Reduce('+', purrr::map(X_iso, ~.x^2)) / length(X_iso),
#       main = "Variance of isotropic process",
#       xlab = "t1", ylab = "t2")
# dev.off()
#

