
bw_smooth <- function(H1, H2, L1, L2, sigma, k) {



}

#' Bi-variate kernel smoothing for functional data on common design
#'
#' @param Y_list List, containing the following elements:
#' - **$t** Vector of sampling points.
#' - **$X** Matrix of observed points, measured on the bi-dimensional grid containing
#' cartesian product of `$t` with itself.
#' @param bw_vec Vector of bandwidths, corresponding to the diagonal elements of
#' the bandwidth matrix.
#' @param xout Vector of evaluation points along one dimension. Smoothing is performed
#' on the cartesian product of `xout` with itself.
#' @returns List, corresponding to the smooth curves.
#' @export


ksmooth_bi <- function(Y_list, bw_vec, xout) {

  tobs <- expand.grid(t1 = Y_list[[1]]$t, t2 = Y_list[[1]]$t)

  tout <- expand.grid(t1 = xout, t2 = xout)

  # Faster than outer
  weights <- vapply(seq_along(tobs$t1),
                    function(x) epa_kernel((tobs$t1[x] - tout$t1) / bw_vec[1])
                    * epa_kernel((tobs$t2[x] - tout$t2) / bw_vec[2]),
                    FUN.VALUE = numeric(length(tout$t1))
                    ) |>
    (\(x) t(sweep(x, MARGIN = 1, STATS = rowSums(x), FUN = "/")))() |>
    (\(x) replace(x, is.nan(x), 0))()

  purrr::map(Y_list,
             ~matrix(crossprod(c(.x$X), weights),
                     nrow = length(xout),
                     ncol = length(xout)
                     )
             )

}



