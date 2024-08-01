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
rout <- 500
alpha_set <- c(pi/30, pi/5, pi/4, pi/3, pi/2 - pi/30)


sigma_set <- c(0.1, 0.5, 1, 2)


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

tic()
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
                                         xout = xout,
                                         sigma = alpha_sheet_sum$sigma_hat)

      # Perform correction for remainder term of g
      alpha_hat_adj <- angle_correct(X_list = sheets_list_sum,
                                     xout = xout,
                                     g_hat = alpha_sheet_sum$g_hat,
                                     alpha_hat = alpha_unique_sum$alpha,
                                     delta = delta,
                                     Hmax = alpha_unique_sum$H_max,
                                     Hmin = alpha_sheet_sum$H_min,
                                     sigma = alpha_sheet_sum$sigma_hat)

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
          'alpha' = alpha_true,
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

toc()

save(result_list,
     file = paste0(here(), "/result_alpha/result_list.rds"))

# analysis of results =========================================================
library(dplyr)
library(tidyr)
library(stringr)
library(pracma)

list_names <- purrr::map_chr(seq_len(nrow(param_cart)),
                             ~paste0("M",  param_cart[.x, "M"],
                                     "_N", param_cart[.x, "N"],
                                    "_alpha", round(param_cart[.x, "alpha"], 2),
                                     "_sigma", param_cart[.x, "sigma"])
)


names(result_list) <- list_names

result_df <- bind_rows(
  lapply(seq_along(result_list), function(result_id) {
    df_result <- as.data.frame(do.call(rbind, result_list[[result_id]]))
    as.data.frame(
      cbind(M = as.numeric(str_extract(list_names[result_id], "(?<=M)[0-9.]+")),
            N = as.numeric(str_extract(list_names[result_id], "(?<=N)[0-9.]+")),
            alpha = as.numeric(str_extract(list_names[result_id],
                                           "(?<=alpha)[0-9.]+")),
            sigma = as.numeric(str_extract(list_names[result_id],
                                           "(?<=sigma)[0-9.]+")),
            alpha_hat = unname(unlist(df_result$alpha_hat)),
            alpha_hat_adj = unname(unlist(df_result$alpha_hat_adj)),
            ghat = unlist(df_result$ghat),
            ghat_adj = unlist(df_result$ghat_adj),
            g_true = unlist(df_result$g_true),
            f_alpha = unname(unlist(df_result$f_alpha))
            )
    )
    }
    )
)


result_long <- result_df |>
  mutate(risk_alpha = abs(alpha - alpha_hat),
         risk_alpha_adj = abs(alpha - alpha_hat_adj),
         risk_ghat = abs(ghat - cot(alpha)),
         risk_ghat_adj = abs(ghat_adj - cot(alpha))) |>
  mutate(risk_rel_adj = risk_alpha_adj / risk_alpha)


result_long |>
  filter(M == 51, N == 100) |>
  select(alpha, sigma, risk_alpha) |>
  ggplot(aes(x = as.factor(sigma)), y = risk_alpha, fill = as.factor(alpha)) +
  geom_boxplot()

lim_plot <- filter(result_df,
                   N == 100,
                   M == 51) |>
  pull(risk_rel) |>
  quantile(c(0, .95))


plot_df(result_df = result_long,
        N = 100,
        M = 51,
        alpha_rem = 0.1,
        var = risk_alpha,
        yaxis = "alpha")

plot_df(result_df = result_long,
        N = 100,
        M = 51,
        alpha_rem = 0.1,
        var = risk_alpha_adj,
        yaxis = "alpha")

 plot_df(result_df = result_long,
        N = 150,
        M = 101,
        alpha_rem = 0.63,
        var = risk_rel_adj,
        yaxis = "alpha") +
   geom_hline(yintercept = 1, col = "red")


plot_df <- function(result_df, N, M, var, yaxis, alpha_rem = NULL) {

  rem_df <- result_df |> filter(alpha %in% alpha_rem)
  var_parse <- rlang::enquo(var)

  df_plot <-  setdiff(x = result_df, y = rem_df) |>
    filter(N == !!N, M == !!M)

  lim_plot <- df_plot |>
    pull(!!var_parse) |>
    quantile(c(0, .99))


  df_plot |>
    ggplot(aes(x = as.factor(sigma), y = !!var_parse,
               fill = as.factor(alpha))) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_continuous(limits = lim_plot) +
    ggtitle(paste0("N = ", N, ", M = ", M, "")) +
    xlab("$\\sigma$") +
    ylab(paste0("$\\mathcal{R}_{", yaxis, "}$")) +
    labs(fill = "$\\alpha$") +
    theme_minimal() +
    scale_fill_grey(start = 0.3, end = 1) +
    theme(plot.title = element_text(hjust = 0.5))

}




