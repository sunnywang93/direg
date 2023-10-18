
# Test to see 1D fbm function works
n_test = 200
end_point = 0.7
plot(seq(0, end_point, length.out = n_test), 
     fbm_fast(0.5, n_test, end_point), type = "l")

test <- lapply(seq_len(20000), function(i) fbm_fft(0.8, 200, 1))

(Reduce('+', lapply(test, function(x) tcrossprod(x))) / length(test)) |>
  (\(x) image(x = seq(0, 1, length.out = 200),
              y = seq(0, 1, length.out = 200),
              z = x))()

par(mfrow = c(1, 2))

H = 0.8
grid_cov = expand.grid(s = seq(0, 1, length.out = n_test), 
                       t = seq(0, 1, length.out = n_test))
bm_cov = 0.5 * (abs(grid_cov$t)**(2*H) + abs(grid_cov$s)**(2*H) - 
                  abs(grid_cov$t - grid_cov$s)**(2*H))
bm_mat = matrix(bm_cov, nrow = n_test, ncol = n_test)
image(bm_mat)

## =====================================================
prod = lapply(seq_len(1000), function(x) fbm_prod(0.8, 0.5, 11, 1))
var_prod <- Reduce('+', lapply(prod, function(x) x^2)) / length(prod)
image(var_prod)


GA::persp3D(seq(0, 1, length.out = 11),
            seq(0, 1, length.out = 11),
            prod[[25]])


prod_rotate <- lapply(seq_len(1000), function(x) 
  fbm_sheet(t1_n = 11, t2_n = 11, e1_n = 11, e2_n = 11,
            H1 = 0.8, H2 = 0.5))

var_rotate <- Reduce('+', lapply(prod_rotate, function(x) x^2)) / 
  length(prod_rotate)

GA::persp3D(seq(0, 1, l = 11),
            seq(0, 1, l = 11),
            var_rotate, theta = -pi/4)

par(mfrow = c(1, 2))
image(var_prod)
image(var_rotate)


## ======================================================
sheets_list <- lapply(seq_len(100), 
                      function(x) {
                        list(
                          t = seq(0, 1, length.out = 20),
                          X = fbm_sheet(20, 20, 20, 20, 0.42, 0.4) 
                        )
                        })


#delta = seq(0.1, 0.5, length.out = 20)
delta = 0.2

alpha_test = estimate_angle(X_list = sheets_list,
                            xout = seq(0, 1, length.out = 11),
                            delta = delta)


alpha_test <- purrr::map(delta,
                         ~estimate_angle(X_list = sheets_list,
                                         xout = seq(0, 1, length.out = 18),
                                         delta = .x))





