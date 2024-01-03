library(direg)
library(foreach)
library(parallel)
library(tictoc)
library(snow)
library(ggplot2)
library(here)

# Parameter settings
Nset <- c(100, 200)
Mset <- c(26, 51)
H1 <- 0.8
H2 <- 0.5
delta_grid <- seq(0.05, 0.4, length.out = 15)
rout <- 500
xout <- seq(0, 1, length.out = 51)
alpha_set <- c(pi/3, pi/4, pi/6, 5*pi/6, pi+pi/3)
delta_c <- 0.25
type_set <- c("sum", "prod")

param_cart <- expand.grid(N = Nset,
                          M = Mset,
                          alpha = alpha_set,
                          type = type_set)
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
                           alpha = pi / 6,
                           H1 = H1,
                           H2 = H2,
                           type = type)
    )

    # Check the variance of the simulated process
    #image(Reduce('+', purrr::map(X_list, ~.x^2)) / length(X_list))

    sheets_list <- purrr::map(X_list,
                              ~list(t = xout,
                                    X = .x))

    alpha_sheet <- estimate_angle(X_list = sheets_list,
                                  xout = xout,
                                  delta = (1 / sqrt(M)) * (1 + delta_c)
                                  )

    alpha_unique <- identify_angle(angles = alpha_sheet,
                                   dout = delta_grid,
                                   xout = xparam)

    id_indic <- names(which.min(abs(alpha_unique - alpha_sheet))) == names(alpha_unique)

    risk <- abs(alpha_unique - alpha_true)

    result <- list(
      'alpha_hat' = alpha_unique,
      'error' = risk,
      'id_indic' = as.double(id_indic)
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


saveRDS(purrr::map_dbl(result_list, ~.x),
        file = paste0("N", N,
                      "_M", M,
                      "alpha", alpha_true,
                      "type", type,
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



