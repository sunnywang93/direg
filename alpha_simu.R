library(direg)
library(foreach)
library(parallel)
library(tictoc)
library(snow)
library(ggplot2)
library(here)

# Parameter settings
N <- 100
M <- 51
H1 <- 0.8
H2 <- 0.5
delta <- 0.2
rout <- 20

xtrue <- seq(0, 1, length.out = M)
xparam <- seq(0, 1, length.out = 21)


# Set seeds to ensure reproducibility
set.seed(123)

seeds <- sample.int(10000,
                    size = rout)

inter_dir <- here('intermediate')

if(!dir.exists(inter_dir)) {
  dir.create(inter_dir)
}
# Set number of cores to use - all available except 2
n_cores <- detectCores() - 2
# Create cluster nodes
cl <- makeCluster(spec = n_cores)
# Register cluster
doSNOW::registerDoSNOW(cl)
# Create progress bar
pb <- txtProgressBar(min = 1,
                     max = rout,
                     style = 3)

opts <- list(progress = function(n) setTxtProgressBar(pb, n))

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
    # Generate fractional brownian sheets
    X_list <- purrr::map(seq_len(N),
                         ~fbm_sheet(
                           t_n = M,
                           e_n = M,
                           alpha_fun = function(x) pi / 6,
                           H1 = H1,
                           H2 = H2)
    )

    # X_can_list <- purrr::map(seq_len(N),
    #                          ~fbm_prod(H1 = H1,
    #                                    H2 = H2,
    #                                    n = M,
    #                                    endpoint = 1)
    # )
    #
    # sheets_can_list <- purrr::map(X_can_list,
    #                           ~list(t = xtrue,
    #                                 X = .x))

    sheets_list <- purrr::map(X_list,
                              ~list(t = xtrue,
                                    X = .x))

    alpha_sheet <- estimate_angle(X_list = sheets_list,
                                  xout = xparam,
                                  delta = 0.2)

    v1 <- c(cos(mean(alpha_sheet$alpha_cot)), sin(mean(alpha_sheet$alpha_cot)))

    v2 <- c(cos(mean(alpha_sheet$alpha_tan)), sin(mean(alpha_sheet$alpha_tan)))

    H_test1 <- H_sheets(X_list = sheets_list,
                        tout = tout,
                        delta = 0.2,
                        base_list = list(v1, v2))


    # Estimate angle between basis vectors
    # alpha_list <- lapply(seq(0.01, 0.5, l = 20), function(x) {
    #   estimate_angle(X_list = sheets_list,
    #                  xout = xparam,
    #                  delta = x)
    # })


    # alpha_mu <- purrr::map_dbl(alpha_list,
    #                        ~mean(.x, na.rm = TRUE)
    #                        )

    save(alpha_mu,
         file = each_filepath)

  }

stopCluster(cl)
toc()

fls <- list.files(inter_dir,
                  pattern = ".rda")


result_list <- lapply(fls,
                      function(x) get(eval(load(paste0(inter_dir, '/', x
                      )))))


saveRDS(purrr::map_dbl(result_list, ~.x),
        file = paste0("N", N,
                      "_M", M,
                      "_H1", H1,
                      "_H2", H2,
                      "delta", delta,
                      ".rds")

)

fs::file_delete(inter_dir)


N100_delta0.1 <- readRDS("N100_M51_H10.8_H20.5delta0.1.rds")
N200_delta0.1 <- readRDS("N200_M51_H10.8_H20.5delta0.1.rds")
N100_delta0.2 <- readRDS("N100_M51_H10.8_H20.5delta0.2.rds")
N200_delta0.2 <- readRDS("N200_M51_H10.8_H20.5delta0.2.rds")

ggplot() +
  geom_boxplot(aes(x = "N100_delta0.1",
                   y = abs(N100_delta0.1 - pi/4),
                   fill = "orange")) +
  geom_boxplot(aes(x = "N200_delta0.1",
                   y = abs(N200_delta0.1 - pi/4),
                   fill = "green")) +
  geom_boxplot(aes(x = "N100_delta0.2",
                   y = abs(N100_delta0.2 - pi/4),
                   fill = "blue")) +
  geom_boxplot(aes(x = "N200_delta0.2",
                   y = abs(N200_delta0.2 - pi/4),
                   fill = "red")) +
  xlab("") +
  ylab("Absolute Errors") +
  ggtitle(paste0(
    "Boxplots for estimated alpha (M =", M,
    ", H1 = ", H1,
    ", H2 = ", H2,
    ", R = ", rout,
    ", True Alpha = ", "pi/4)")
  ) +
  theme(legend.position = "none")



