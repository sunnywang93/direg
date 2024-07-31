library(direg)
library(foreach)
library(parallel)
library(tictoc)
library(snow)
library(ggplot2)
library(here)

# Parameter settings
Nset <- c(100, 150)
Mset <- c(51, 101)
H1 <- 0.8
H2 <- 0.5
rout <- 50
alpha_set <- c(pi/30, pi/5, pi/4, pi/3, pi/2 - pi/30)


sigma_set <- c(0.05, 0.1)


param_cart <- expand.grid(M = Mset,
                          N = Nset,
                          alpha = alpha_set,
                          sigma = sigma_set)
# Set seeds to ensure reproducibility
set.seed(123)

seeds <- sample.int(10000,
                    size = rout)

result_folder <- here('result_alpha')
if(!dir.exists(result_folder)) {
  dir.create(result_folder)
}

# Set number of cores to use
n_cores <- 50
# Create cluster nodes
cl <- makeCluster(spec = n_cores)
# Register cluster
doSNOW::registerDoSNOW(cl)

result_list <- foreach(k = seq_len(nrow(param_cart))) %do% {

  alpha_true <- param_cart[k, "alpha"]
  M <- param_cart[k, "M"]
  N <- param_cart[k, "N"]
  delta_grid <- seq(1/sqrt(M), 0.4, length.out = 15)
  sigma <- param_cart[k, "sigma"]

  foreach(i = 1:rout,
          .packages = c("direg"),
          .combine = 'c',
          .multicombine = TRUE) %dopar%
    {

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

      list(
        list(
          'alpha_hat' = alpha_unique_sum$alpha,
          'alpha_hat_adj' = alpha_hat_adj$alpha_adj,
          'ghat' = alpha_sheet_sum$g_hat,
          'ghat_adj' = alpha_hat_adj$g_adj,
          'g_true' = g_true,
          'f_alpha' = alpha_hat_adj$f_alpha
        )
      )
    }

}

save(result_list,
     file = paste0(here(), "/result_alpha/result_list.rds"))

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

# need to compute against cotangent of alpha, not g_true!
result_df <- result_df |>
  mutate(risk_alpha_adj = abs(alpha - alpha_hat_adj),
         risk_alpha_hat = abs(alpha - alpha_hat),
         risk_ghat = abs(ghat - cot(alpha)),
         risk_ghat_adj = abs(ghat_adj - cot(alpha))) |>
  mutate(risk_rel = risk_alpha_hat / risk_alpha_adj)

# result_long <- result_long |>
#   mutate(risk = risk - pi/2)
#
#
# result_long <- result_long |>
#   mutate(risk = abs(risk))

plot_set <- expand.grid(N = Nset, M = Mset)

alpha_adj_list <- lapply(seq_len(nrow(plot_set)), function(row) {
  plot_df(result_df = result_df,
          N = plot_set[row, "N"],
          M = plot_set[row, "M"],
          alpha_rem = c(2.09, 2.36, 2.62),
          var = risk_alpha_adj,
          yaxis = "alpha")
})


alpha_adj <- wrap_plots(alpha_adj_list, nrow = 2, ncol = 2, guides = "collect")

for(i in 1:length(alpha_adj_list)) {

  tikz(file = paste0("/home/swang/directional_regularity/plots/alpha_adj/",
                     alpha_adj_list[[i]]$labels$title, ".tex"),
       width = 5,
       height = 5,
       standAlone = TRUE)

  plot(alpha_adj_list[[i]])

  dev.off()

}



alpha_list <- lapply(seq_len(nrow(plot_set)), function(row) {
  plot_df(result_df = result_df,
          N = plot_set[row, "N"],
          M = plot_set[row, "M"],
          alpha_rem = c(2.09, 2.36, 2.62),
          var = risk_alpha_hat,
          yaxis = "alpha")
})

alpha_hat <- wrap_plots(alpha_list, nrow = 2, ncol = 2, guides = "collect")


for(i in 1:length(alpha_list)) {

  tikz(file = paste0("/home/swang/directional_regularity/plots/alpha/",
                     alpha_list[[i]]$labels$title, ".tex"),
       width = 5,
       height = 5,
       standAlone = TRUE)

  plot(alpha_list[[i]])

  dev.off()

}





alpha_rel_list <- lapply(seq_len(nrow(plot_set)), function(row) {

  lim_plot <- filter(result_df,
                     N == plot_set[row, "N"],
                     M == plot_set[row, "M"]) |>
    pull(risk_rel) |>
    quantile(c(0, .95))

  plot_df(result_df = result_df,
          N = plot_set[row, "N"],
          M = plot_set[row, "M"],
          alpha_rem = c(2.09, 2.36, 2.62),
          var = risk_rel,
          yaxis = "alpha") +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(limits = lim_plot) +
    geom_hline(yintercept = 1, col = "red")
})

for(i in 1:length(alpha_rel_list)) {

  tikz(file = paste0("/home/swang/directional_regularity/plots/alpha_rel/",
                     alpha_rel_list[[i]]$labels$title, ".tex"),
       width = 5,
       height = 5,
       standAlone = TRUE)

  plot(alpha_rel_list[[i]])

  dev.off()

}



alpha_rel <- wrap_plots(alpha_rel_list, nrow = 2, ncol = 2, guides = "collect")


f_list <- lapply(seq_len(nrow(plot_set)), function(row) {
  plot_df(result_df = result_df,
          N = plot_set[row, "N"],
          M = plot_set[row, "M"],
          alpha_rem = c(2.09, 2.36, 2.62),
          var = f_alpha,
          yaxis = "f")
})

f_rel <- wrap_plots(f_list, nrow = 2, ncol = 2, guides = "collect")


ghat_adj_list <- lapply(seq_len(nrow(plot_set)), function(row) {
  plot_df(result_df = result_df,
          N = plot_set[row, "N"],
          M = plot_set[row, "M"],
          alpha_rem = c(0.1, 2.09, 2.36, 2.62),
          var = risk_ghat_adj,
          yaxis = "g_adj")
})

ghat_adj <- wrap_plots(ghat_adj_list, nrow = 2, ncol = 2, guides = "collect")


for(i in 1:length(ghat_adj_list)) {

  tikz(file = paste0("/home/swang/directional_regularity/plots/g_adj/",
                     ghat_adj_list[[i]]$labels$title, ".tex"),
       width = 5,
       height = 5,
       standAlone = TRUE)

  plot(ghat_adj_list[[i]])

  dev.off()

}


ghat_list <- lapply(seq_len(nrow(plot_set)), function(row) {
  plot_df(result_df = result_df,
          N = plot_set[row, "N"],
          M = plot_set[row, "M"],
          alpha_rem = c(0.1, 2.09, 2.36, 2.62),
          var = risk_ghat,
          yaxis = "g")
})

ghat <- wrap_plots(ghat_list, nrow = 2, ncol = 2, guides = "collect")


for(i in 1:length(ghat_list)) {

  tikz(file = paste0("/home/swang/directional_regularity/plots/g/",
                     ghat_list[[i]]$labels$title, ".tex"),
       width = 5,
       height = 5,
       standAlone = TRUE)

  plot(ghat_list[[i]])

  dev.off()

}

# lim_N100_M51_prod <- filter(result_long,
#                             dgp == "prod",
#                             N == 100,
#                             M == 51) |>
#   pull(risk) |>
#   quantile(c(.1, .9))
#
# boxplot_N100_M51_prod <- result_long |>
#   ggplot(aes(x = as.factor(sigma), y = risk,
#              fill = as.factor(alpha))) +
#   geom_boxplot(outlier.shape = NA) +
#   scale_y_continuous(limits = quantile(lim_N100_M51_prod, c(0.1, 0.9))) +
#   ggtitle("Boxplots $f_1 = B_1 * B2 (N = 100, M = 51)$") +
#   xlab("$\\sigma$") +
#   ylab("$\\mathcal{R}_{\\alpha}$") +
#   labs(fill = "$\\alpha$") +
#   theme_minimal() +
#   scale_fill_grey(start = 0.3, end = 1)


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


plot_df <- function(result_df, N, M, var, yaxis, alpha_rem = NULL) {

  rem_df <- result_df |> filter(alpha %in% alpha_rem)
  var_parse <- rlang::enquo(var)

  setdiff(x = result_df, y = rem_df) |>
    filter(N == !!N, M == !!M) |>
    ggplot(aes(x = as.factor(sigma), y = !!var_parse,
               fill = as.factor(alpha))) +
    geom_boxplot() +
    ggtitle(paste0("N = ", N, ", M = ", M, "")) +
    xlab("$\\sigma$") +
    ylab(paste0("$\\mathcal{R}_{", yaxis, "}$")) +
    labs(fill = "$\\alpha$") +
    theme_minimal() +
    scale_fill_grey(start = 0.3, end = 1) +
    theme(plot.title = element_text(hjust = 0.5))

}





